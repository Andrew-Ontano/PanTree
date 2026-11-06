import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import pysam
import pantreelib
from scripts.calculate import generate_window_tasks


def run(args):
    logging.info(f"Preparing window tasks for local module from input VCF: {args.input_vcf}")

    # Inspect VCF header for sample list
    vcf = pysam.VariantFile(args.input_vcf)
    all_vcf_samples = list(vcf.header.samples)
    target_samples = args.samples.split(",") if args.samples else None

    if target_samples:
        samples = [s for s in all_vcf_samples if s in target_samples]
    else:
        samples = all_vcf_samples

    active_samples = ["REF"] + samples if args.include_ref else samples
    vcf.close()

    tasks = generate_window_tasks(args)
    logging.info(f"Generated {len(tasks)} window task(s). Processing with {args.threads} thread(s)...")

    # Output columns: Chromosome, Window_start, Window_end, Total_variants, Called_variants, Group, Count
    header = [
        "Chromosome", "Window_start", "Window_end",
        "Total_variants", "Called_variants", "Group", "Count"
    ]

    all_rows = []

    if args.threads > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            future_to_task = {
                executor.submit(pantreelib.process_local_window_task, task): task for task in tasks
            }
            count = 0
            for future in as_completed(future_to_task):
                count += 1
                if args.verbose and count % 10 == 0:
                    logging.debug(f"Processed {count}/{len(tasks)} windows...")
                try:
                    res_rows = future.result()
                    all_rows.extend(res_rows)
                except Exception as exc:
                    logging.error(f"Task generated an exception: {exc}")
    else:
        for idx, task in enumerate(tasks, start=1):
            if args.verbose and idx % 10 == 0:
                logging.debug(f"Processing window {idx}/{len(tasks)}...")
            res_rows = pantreelib.process_local_window_task(task)
            all_rows.extend(res_rows)

    # Mapping samples to numeric IDs: 0 for REF, 1 for sample1, ..., N for sampleN
    # If REF is not included, sample1..sampleN are numbered 1..N
    sample_to_id = {}
    id_to_sample = {}
    if args.include_ref:
        sample_to_id["REF"] = 0
        id_to_sample[0] = "REF"
        for idx, s in enumerate(samples, start=1):
            sample_to_id[s] = idx
            id_to_sample[idx] = s
    else:
        for idx, s in enumerate(samples, start=1):
            sample_to_id[s] = idx
            id_to_sample[idx] = s

    # Convert Group sample names to numeric IDs in output rows
    for row in all_rows:
        sample_names = row["Group"].split(",")
        id_list = [str(sample_to_id[s]) for s in sample_names if s in sample_to_id]
        row["Group_numeric"] = ",".join(id_list)

    # Sorting order: Chromosome > Window_start > Group (by numeric group / VCF order)
    all_rows.sort(key=lambda r: (r["Chromosome"], r["Window_start"], [int(x) for x in r["Group_numeric"].split(",")]))

    # Write output TSV with commented section above the table
    logging.info(f"Writing results to {args.output_tsv}...")
    with open(args.output_tsv, "w") as out_f:
        out_f.write("# PanTree local: Sample ID Mapping\n")
        for num_id in sorted(id_to_sample.keys()):
            out_f.write(f"# ID {num_id}: {id_to_sample[num_id]}\n")
        out_f.write("\t".join(header) + "\n")
        for row in all_rows:
            row_dict = dict(row)
            row_dict["Group"] = row["Group_numeric"]
            row_vals = [str(row_dict.get(col, "NA")) for col in header]
            out_f.write("\t".join(row_vals) + "\n")

    logging.info("Local module execution complete!")
