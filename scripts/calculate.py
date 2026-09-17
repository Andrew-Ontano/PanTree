import logging
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
import pysam
import pantreelib


def generate_window_tasks(args):
    vcf = pysam.VariantFile(args.input_vcf)
    target_samples = args.samples.split(",") if args.samples else None

    # Retrieve contig lengths and names
    chrom_lens = {contig: record.length for contig, record in vcf.header.contigs.items()}

    step = args.window_size - args.overlap if args.overlap > 0 else args.window_size
    step = max(1, step)

    tasks = []

    if args.fixed_variants:
        for chrom in vcf.header.contigs:
            try:
                records = vcf.fetch(contig=chrom)
            except ValueError:
                continue
            pos_list = [r.pos for r in records]
            if not pos_list:
                continue
            num_vars = len(pos_list)
            for i in range(0, num_vars, step):
                start = pos_list[i]
                end_idx = min(i + args.window_size - 1, num_vars - 1)
                end = pos_list[end_idx]
                tasks.append((args.input_vcf, chrom, start, end, target_samples, args.include_ref))
    else:
        for chrom, length in chrom_lens.items():
            if not length or length <= 0:
                # If length not in header, try fetching all records for contig to find max pos
                try:
                    records = list(vcf.fetch(contig=chrom))
                    if records:
                        length = records[-1].pos
                    else:
                        continue
                except ValueError:
                    continue

            for start in range(1, length + 1, step):
                end = min(start + args.window_size - 1, length)
                tasks.append((args.input_vcf, chrom, start, end, target_samples, args.include_ref))

    vcf.close()
    return tasks


def run(args):
    logging.info(f"Preparing window tasks for input VCF: {args.input_vcf}")

    # Inspect VCF header for sample list and columns
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

    # Define header columns
    header = [
        "Chromosome", "Window_start", "Window_end",
        "Total_variants", "Total_SNPs", "Total_MNPs", "Total_INDELs",
        "Total_bp", "Tree_newick"
    ]

    for s in active_samples:
        header.extend([
            f"{s}_called_sites", f"{s}_snps", f"{s}_mnps", f"{s}_indels",
            f"{s}_bp", f"{s}_unique_alleles"
        ])

    pairs = list(itertools.combinations(active_samples, 2))
    for s1, s2 in pairs:
        header.extend([
            f"{s1}-{s2}_patristic", f"{s1}-{s2}_jaccard", f"{s1}-{s2}_manhattan"
        ])

    results = []

    if args.threads > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            future_to_task = {
                executor.submit(pantreelib.process_window_task, task): task for task in tasks
            }
            count = 0
            for future in as_completed(future_to_task):
                count += 1
                if args.verbose and count % 10 == 0:
                    logging.debug(f"Processed {count}/{len(tasks)} windows...")
                try:
                    res = future.result()
                    results.append(res)
                except Exception as exc:
                    logging.error(f"Task generated an exception: {exc}")
    else:
        for idx, task in enumerate(tasks, start=1):
            if args.verbose and idx % 10 == 0:
                logging.debug(f"Processing window {idx}/{len(tasks)}...")
            res = pantreelib.process_window_task(task)
            results.append(res)

    # Sort results by Chromosome and Window_start
    results.sort(key=lambda r: (r["Chromosome"], r["Window_start"]))

    # Write output TSV
    logging.info(f"Writing results to {args.output_tsv}...")
    with open(args.output_tsv, "w") as out_f:
        out_f.write("\t".join(header) + "\n")
        for res in results:
            row_vals = [str(res.get(col, "NA")) for col in header]
            out_f.write("\t".join(row_vals) + "\n")

    logging.info("Calculation complete!")
