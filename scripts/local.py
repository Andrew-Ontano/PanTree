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

    # Sorting order: Chromosome > Window_start > Group
    all_rows.sort(key=lambda r: (r["Chromosome"], r["Window_start"], r["Group"]))

    # Write output TSV
    logging.info(f"Writing results to {args.output_tsv}...")
    with open(args.output_tsv, "w") as out_f:
        out_f.write("\t".join(header) + "\n")
        for row in all_rows:
            row_vals = [str(row.get(col, "NA")) for col in header]
            out_f.write("\t".join(row_vals) + "\n")

    logging.info("Local module execution complete!")
