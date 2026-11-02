import io
import itertools
import logging
import numpy as np
import pysam
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def parse_sample_genotypes(variant, sample_list):
    """
    Parse sample genotypes from a pysam VariantRecord.
    Returns:
      gts: dict mapping sample_name -> allele_index (int) or None if missing/uncalled.
    """
    gts = {}
    for s in sample_list:
        gt = variant.samples[s].get('GT')
        if gt is None or None in gt:
            gts[s] = None
        else:
            non_ref = [a for a in gt if a > 0]
            if non_ref:
                gts[s] = non_ref[0]
            else:
                gts[s] = gt[0]
    return gts


def compute_sequence_distance(seq1, seq2):
    """
    Compute p-distance (proportion of nucleotide differences) between two aligned sequences.
    Ignores positions where either sequence has a missing/gap '-'.
    """
    valid_len = 0
    diffs = 0
    for char1, char2 in zip(seq1, seq2):
        if char1 != '-' and char2 != '-':
            valid_len += 1
            if char1 != char2:
                diffs += 1

    if valid_len == 0:
        return 0.0
    return diffs / valid_len


def process_window_task(args):
    """
    Worker task to process a single window.
    args: tuple (vcf_path, chrom, start, end, target_samples, include_ref)
    """
    vcf_path, chrom, start, end, target_samples, include_ref = args

    vcf = pysam.VariantFile(vcf_path)
    all_vcf_samples = list(vcf.header.samples)

    if target_samples:
        samples = [s for s in all_vcf_samples if s in target_samples]
    else:
        samples = all_vcf_samples

    active_samples = ["REF"] + samples if include_ref else samples

    result = {
        "Chromosome": chrom,
        "Window_start": start,
        "Window_end": end,
        "Total_variants": 0,
        "Total_SNPs": 0,
        "Total_MNPs": 0,
        "Total_INDELs": 0,
        "Total_bp": 0,
        "Tree_newick": "NA"
    }

    # Initialize sample-specific fields
    for s in active_samples:
        result[f"{s}_called_sites"] = 0
        result[f"{s}_snps"] = 0
        result[f"{s}_mnps"] = 0
        result[f"{s}_indels"] = 0
        result[f"{s}_bp"] = 0
        result[f"{s}_unique_alleles"] = 0

    pairs = list(itertools.combinations(active_samples, 2))
    for s1, s2 in pairs:
        result[f"{s1}-{s2}_patristic"] = "NA"
        result[f"{s1}-{s2}_jaccard"] = "NA"
        result[f"{s1}-{s2}_manhattan"] = "NA"

    # Fetch variant records
    try:
        records = list(vcf.fetch(contig=chrom, start=start - 1, stop=end))
    except ValueError:
        records = []
    vcf.close()

    if not records:
        return result

    # Store sequence alignment characters per sample: list of characters across SNP sites
    sample_seq_chars = {s: [] for s in active_samples}

    # Matrices for genotype distance calculations (sample x variant)
    gt_matrix = []
    alt_bool_matrix = []

    for variant in records:
        ref = variant.ref
        alts = list(variant.alts or [])
        if not alts:
            continue

        # Expect supplied variants are SNPs (len(ref) == 1 and all len(alts) == 1)
        if len(ref) != 1 or any(len(a) != 1 for a in alts):
            continue

        result["Total_variants"] += 1
        result["Total_SNPs"] += 1

        allele_dict = {0: ref}
        for idx, alt in enumerate(alts, start=1):
            allele_dict[idx] = alt

        gts = parse_sample_genotypes(variant, samples)
        if include_ref:
            gts["REF"] = 0

        # Unique allele check across active_samples
        called_alleles = [gts[s] for s in active_samples if gts[s] is not None]
        allele_counts = {}
        for a in called_alleles:
            allele_counts[a] = allele_counts.get(a, 0) + 1

        var_gt_row = []
        var_alt_row = []

        for s in active_samples:
            allele_idx = gts[s]
            if allele_idx is not None:
                result[f"{s}_called_sites"] += 1
                result[f"{s}_snps"] += 1

                if allele_counts.get(allele_idx, 0) == 1:
                    result[f"{s}_unique_alleles"] += 1

                nuc = allele_dict.get(allele_idx, '-')
                sample_seq_chars[s].append(nuc)
                result[f"{s}_bp"] += 1

                var_gt_row.append(allele_idx)
                var_alt_row.append(1 if allele_idx > 0 else 0)
            else:
                sample_seq_chars[s].append('-')
                var_gt_row.append(np.nan)
                var_alt_row.append(np.nan)

        gt_matrix.append(var_gt_row)
        alt_bool_matrix.append(var_alt_row)

    if result["Total_variants"] == 0:
        return result

    # Total window sequence length equals number of SNP sites
    result["Total_bp"] = result["Total_variants"]
    full_sample_seqs = {s: "".join(sample_seq_chars[s]) for s in active_samples}

    # --- 1. Compute Pairwise Jaccard and Manhattan Distances ---
    gt_arr = np.array(gt_matrix).T
    alt_arr = np.array(alt_bool_matrix).T
    num_samples = len(active_samples)

    manhattan_dist = np.zeros((num_samples, num_samples))
    jaccard_dist = np.full((num_samples, num_samples), np.nan)

    for i in range(num_samples):
        for j in range(i + 1, num_samples):
            s1_name = active_samples[i]
            s2_name = active_samples[j]

            # Manhattan distance normalized by valid (non-missing) sites
            valid_mask = ~np.isnan(gt_arr[i]) & ~np.isnan(gt_arr[j])
            valid_count = np.sum(valid_mask)
            if valid_count > 0:
                diff = np.sum(np.abs(gt_arr[i, valid_mask] - gt_arr[j, valid_mask]))
                norm_diff = diff / valid_count
                manhattan_dist[i, j] = norm_diff
                manhattan_dist[j, i] = norm_diff
                result[f"{s1_name}-{s2_name}_manhattan"] = round(norm_diff, 6)
            else:
                result[f"{s1_name}-{s2_name}_manhattan"] = "NA"

            # Jaccard distance on ALT presence (omit REF comparisons if include_ref is True)
            if include_ref and (s1_name == "REF" or s2_name == "REF"):
                result[f"{s1_name}-{s2_name}_jaccard"] = "NA"
            else:
                v1 = alt_arr[i, valid_mask]
                v2 = alt_arr[j, valid_mask]
                if len(v1) > 0:
                    both_alt = np.sum((v1 == 1) & (v2 == 1))
                    either_alt = np.sum((v1 == 1) | (v2 == 1))
                    if either_alt > 0:
                        j_dist = 1.0 - (both_alt / either_alt)
                        jaccard_dist[i, j] = j_dist
                        jaccard_dist[j, i] = j_dist
                        result[f"{s1_name}-{s2_name}_jaccard"] = round(j_dist, 6)
                    else:
                        result[f"{s1_name}-{s2_name}_jaccard"] = 0.0

    # --- 2. Compute Sequence Alignment Distance Matrix & Build NJ Tree ---
    seq_dist_matrix = np.zeros((num_samples, num_samples))
    for i in range(num_samples):
        for j in range(i + 1, num_samples):
            dist = compute_sequence_distance(
                full_sample_seqs[active_samples[i]],
                full_sample_seqs[active_samples[j]]
            )
            seq_dist_matrix[i, j] = dist
            seq_dist_matrix[j, i] = dist

    if num_samples >= 3:
        lower_tri = []
        for i in range(num_samples):
            lower_tri.append(seq_dist_matrix[i, :i + 1].tolist())

        try:
            dm = DistanceMatrix(active_samples, lower_tri)
            constructor = DistanceTreeConstructor()
            nj_tree = constructor.nj(dm)

            handle = io.StringIO()
            Phylo.write(nj_tree, handle, "newick")
            result["Tree_newick"] = handle.getvalue().strip()

            for s1, s2 in pairs:
                try:
                    patristic = nj_tree.distance(s1, s2)
                    result[f"{s1}-{s2}_patristic"] = round(patristic, 6)
                except Exception:
                    result[f"{s1}-{s2}_patristic"] = "NA"
        except Exception:
            result["Tree_newick"] = "TREE_ERROR"
    elif num_samples == 2:
        s1, s2 = active_samples[0], active_samples[1]
        d = seq_dist_matrix[0, 1]
        half_d = round(d / 2.0, 6)
        result["Tree_newick"] = f"({s1}:{half_d},{s2}:{half_d});"
        result[f"{s1}-{s2}_patristic"] = round(d, 6)

    return result


def process_local_window_task(args):
    """
    Worker task to process a single window for the 'local' module.
    args: tuple (vcf_path, chrom, start, end, target_samples, include_ref)

    Returns a list of dicts with keys:
      Chromosome, Window_start, Window_end, Total_variants, Called_variants, Group, Count
    """
    vcf_path, chrom, start, end, target_samples, include_ref = args

    vcf = pysam.VariantFile(vcf_path)
    all_vcf_samples = list(vcf.header.samples)

    if target_samples:
        samples = [s for s in all_vcf_samples if s in target_samples]
    else:
        samples = all_vcf_samples

    active_samples = ["REF"] + samples if include_ref else samples

    # Fetch variant records
    try:
        records = list(vcf.fetch(contig=chrom, start=start - 1, stop=end))
    except ValueError:
        records = []
    vcf.close()

    total_variants = 0
    # Store parsed called genotypes per site: list of dicts {sample: allele_index}
    site_gts_list = []

    for variant in records:
        ref = variant.ref
        alts = list(variant.alts or [])
        if not alts:
            continue

        # Expect supplied variants are SNPs (len(ref) == 1 and all len(alts) == 1)
        if len(ref) != 1 or any(len(a) != 1 for a in alts):
            continue

        total_variants += 1

        gts = parse_sample_genotypes(variant, samples)
        if include_ref:
            gts["REF"] = 0

        site_gts_list.append(gts)

    if total_variants == 0:
        return []

    # Count polarized groupings across all sites in this window
    # group_key is a tuple of sample names in active_samples order
    group_counts = {}

    for gts in site_gts_list:
        # Non-missing samples at this site
        called_samples = [s for s in active_samples if gts.get(s) is not None]
        if not called_samples:
            continue

        # Map allele -> list of samples carrying that allele
        allele_to_samples = {}
        for s in called_samples:
            allele = gts[s]
            allele_to_samples.setdefault(allele, []).append(s)

        # A grouping is formed if all samples in the grouping share an allele
        # that is unique against all other non-missing evaluated samples at that site.
        # This corresponds exactly to each distinct allele present at this site!
        # If all samples share allele A, it forms a group of size N with that allele.
        # If samples partition into alleles A and B, they form two groups.
        for allele, grp in allele_to_samples.items():
            # Group key keeps samples in active_samples order
            grp_sorted = tuple(s for s in active_samples if s in grp)
            group_counts[grp_sorted] = group_counts.get(grp_sorted, 0) + 1

    # If no groupings have Count > 0, return empty list
    if not group_counts:
        return []

    results = []
    for grp_tuple, count in group_counts.items():
        if count <= 0:
            continue

        grp_set = set(grp_tuple)
        # Called_variants: number of sites in the window where all samples in the grouping had a called genotype
        called_variants = 0
        for gts in site_gts_list:
            if all(gts.get(s) is not None for s in grp_set):
                called_variants += 1

        group_str = ",".join(grp_tuple)
        results.append({
            "Chromosome": chrom,
            "Window_start": start,
            "Window_end": end,
            "Total_variants": total_variants,
            "Called_variants": called_variants,
            "Group": group_str,
            "Count": count
        })

    return results
