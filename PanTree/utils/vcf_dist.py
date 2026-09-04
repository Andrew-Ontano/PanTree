#!/usr/bin/env python3
import argparse
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.spatial.distance import pdist, squareform
import pysam
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import io


def parse_pysam_genotypes(variant, sample_list):
    """
    Parses genotypes from pysam VariantRecord into ALT allele counts.
    Handles haploid/diploid data seamlessly (0/0 -> 0, 0/1 -> 1, 1/1 -> 2, 1 -> 1).
    """
    parsed = []
    for sample in sample_list:
        gt = variant.samples[sample].get('GT')
        # If missing data (.)
        if gt is None or None in gt:
            parsed.append(np.nan)
        else:
            parsed.append(sum(gt))
    return np.array(parsed, dtype=np.float32)


def process_window(args):
    """Worker function to process a single genomic window using pysam."""
    vcf_path, chrom, start, end, samples, include_ref = args
    # Open VCF with pysam (each process needs its own file handle)
    vcf = pysam.VariantFile(vcf_path)
    # Get active sample names from header
    vcf_samples = list(vcf.header.samples)
    if samples:
        vcf_samples = [s for s in vcf_samples if s in samples]

    active_samples = ["REF"] + vcf_samples if include_ref else vcf_samples

    # Initialize result dictionary with NAs
    pairs = list(itertools.combinations(active_samples, 2))
    row_result = {
        "Chromosome": chrom,
        "Window_start": start,
        "Window_end": end,
        "Variant_count": 0,
        "Tree_newick": "NA"
    }

    for (s1, s2) in pairs:
        row_result[f"{s1}-{s2}_tree"] = "NA"
        row_result[f"{s1}-{s2}_jaccard"] = "NA"
        row_result[f"{s1}-{s2}_d"] = "NA"

    gt_matrix = []

    # pysam fetch uses 0-based half-open coordinates [start-1, end)
    try:
        records = vcf.fetch(contig=chrom, start=start - 1, stop=end)
    except ValueError:
        # Contig not found in index or empty region
        records = []

    for variant in records:
        # Filter for biallelic SNPs
        if len(variant.alts or []) != 1 or len(variant.ref) != 1 or len(variant.alts[0]) != 1:
            continue

        gt = parse_pysam_genotypes(variant, vcf_samples)
        if include_ref:
            gt = np.insert(gt, 0, 0.0)  # Prepend '0' for the reference

        gt_matrix.append(gt)

    vcf.close()

    var_count = len(gt_matrix)
    row_result["Variant_count"] = var_count

    # If no variants in window, return early with NAs
    if var_count == 0:
        return row_result

    # Convert to array and transpose to (samples x variants)
    gt_matrix = np.array(gt_matrix).T

    # Fill NaNs with 0 for scipy's pdist compatibility (assuming missing == Ref for quick math)
    # Note: If your data has very high missingness, you may want a custom pairwise function here.
    gt_matrix = np.nan_to_num(gt_matrix, nan=0.0)

    # --- 1. Jaccard & Pairwise 'd' Distances ---
    # Convert genotypes > 0 to 1 for boolean jaccard
    bool_matrix = (gt_matrix > 0).astype(int)
    jaccard_dists = pdist(bool_matrix, metric='jaccard')
    jaccard_sq = squareform(jaccard_dists)

    # Absolute Pairwise Distance (d) - Manhattan distance normalized by valid sites
    d_dists = pdist(gt_matrix, metric='cityblock') / var_count
    d_sq = squareform(d_dists)

    # --- 2. Neighbor-Joining Tree & Patristic Distances ---
    lower_tri = []
    for i in range(len(active_samples)):
        lower_tri.append(d_sq[i, :i + 1].tolist())

    dm = DistanceMatrix(active_samples, lower_tri)
    constructor = DistanceTreeConstructor()

    try:
        nj_tree = constructor.nj(dm)

        # Get Newick string
        handle = io.StringIO()
        Phylo.write(nj_tree, handle, "newick")
        row_result["Tree_newick"] = handle.getvalue().strip()

        # --- 3. Compile Pairwise Results ---
        for i, (s1, s2) in enumerate(pairs):
            idx1 = active_samples.index(s1)
            idx2 = active_samples.index(s2)

            # Extract patristic distance directly from the NJ tree
            patristic = nj_tree.distance(s1, s2)

            row_result[f"{s1}-{s2}_tree"] = round(patristic, 6)
            row_result[f"{s1}-{s2}_jaccard"] = round(jaccard_sq[idx1, idx2], 6)
            row_result[f"{s1}-{s2}_d"] = round(d_sq[idx1, idx2], 6)

    except Exception as e:
        row_result["Tree_newick"] = "TREE_ERROR"

    return row_result


def main():
    parser = argparse.ArgumentParser(
        description="Build NJ trees and calculate distances for windowed VCFs (pysam backend).")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file (indexed with tabix or csi)")
    parser.add_argument("-w", "--window", type=int, default=50000, help="Window size in bp")
    parser.add_argument("-s", "--samples", type=str, help="Comma-separated list of samples to include")
    parser.add_argument("-r", "--include-ref", action="store_true", help="Include the reference genome in comparisons")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for parallel processing")
    args = parser.parse_args()

    samples_list = args.samples.split(",") if args.samples else None

    # Read VCF header to get chromosomes and lengths
    vcf = pysam.VariantFile(args.vcf)

    vcf_samples = list(vcf.header.samples)
    if samples_list:
        vcf_samples = [s for s in vcf_samples if s in samples_list]

    active_samples = ["REF"] + vcf_samples if args.include_ref else vcf_samples

    # Extract contig lengths directly from the VCF header
    chrom_lens = {contig: record.length for contig, record in vcf.header.contigs.items()}
    vcf.close()

    print(chrom_lens)

    # Generate window coordinates
    tasks = []
    for chrom, length in chrom_lens.items():
        # If length is None (missing from header), skip or handle accordingly
        if not length:
            continue
        for start in range(1, length + 1, args.window):
            end = min(start + args.window - 1, length)
            tasks.append((args.vcf, chrom, start, end, samples_list, args.include_ref))

    # Prepare TSV Header
    pairs = list(itertools.combinations(active_samples, 2))
    header = ["Chromosome", "Window_start", "Window_end", "Variant_count", "Tree_newick"]
    for s1, s2 in pairs:
        header.extend([f"{s1}-{s2}_tree", f"{s1}-{s2}_jaccard", f"{s1}-{s2}_d"])

    # Process in parallel
    print(f"Processing {len(tasks)} windows across {args.threads} threads...")

    with open(args.output, 'w') as out_f:
        out_f.write("\t".join(header) + "\n")

        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            future_to_window = {executor.submit(process_window, task): task for task in tasks}

            for future in as_completed(future_to_window):
                try:
                    res = future.result()
                    row_vals = [str(res.get(col, "NA")) for col in header]
                    out_f.write("\t".join(row_vals) + "\n")
                except Exception as exc:
                    print(f"Window generated an exception: {exc}")

    print(f"Complete! Results written to {args.output}")


if __name__ == "__main__":
    main()