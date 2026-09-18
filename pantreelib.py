import io
import itertools
import logging
import tempfile
import subprocess
import os
import numpy as np
import pysam
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import SeqIO


def classify_variant(ref, alts):
    """
    Classify a variant as SNP, MNP, or INDEL.
    ref: string
    alts: list of strings (or tuple)
    """
    all_alleles = [ref] + list(alts)
    ref_len = len(ref)
    allele_lens = [len(a) for a in all_alleles]

    if all(l == 1 for l in allele_lens):
        return "SNP"
    elif all(l == ref_len for l in allele_lens) and ref_len > 1:
        return "MNP"
    else:
        return "INDEL"


def align_alleles_with_mafft(allele_dict):
    """
    Given a dict mapping allele_index (0, 1, 2, ...) to sequence string,
    run MAFFT to align them if lengths differ.
    Returns dict mapping allele_index to aligned sequence string (with gaps '-').
    """
    seq_lengths = {len(seq) for seq in allele_dict.values()}
    if len(seq_lengths) == 1:
        # Prealigned, same length
        return allele_dict

    # Construct FASTA input safely with tempfile
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as in_f:
        in_path = in_f.name
        for idx, seq in allele_dict.items():
            in_f.write(f">{idx}\n{seq}\n")

    aligned_dict = {}
    try:
        cmd = ["mafft", "--auto", "--quiet", in_path]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0 and res.stdout:
            handle = io.StringIO(res.stdout)
            for record in SeqIO.parse(handle, "fasta"):
                idx = int(record.id)
                aligned_dict[idx] = str(record.seq).upper()
        else:
            # Fallback if MAFFT fails: pad right with '-'
            max_len = max(seq_lengths)
            for idx, seq in allele_dict.items():
                aligned_dict[idx] = seq.ljust(max_len, '-')
    except Exception:
        # Fallback padding if executable not found
        max_len = max(seq_lengths)
        for idx, seq in allele_dict.items():
            aligned_dict[idx] = seq.ljust(max_len, '-')
    finally:
        if os.path.exists(in_path):
            os.remove(in_path)

    return aligned_dict


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
    Ignores positions where either sequence has a gap '-'.
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

    # Store sequence alignment components per sample: list of sequence strings for each variant
    sample_alignments = {s: [] for s in active_samples}

    # Matrices for genotype distance calculations (sample x variant)
    gt_matrix = []
    alt_bool_matrix = []

    for variant in records:
        ref = variant.ref
        alts = list(variant.alts or [])
        if not alts:
            continue

        var_type = classify_variant(ref, alts)
        result["Total_variants"] += 1
        if var_type == "SNP":
            result["Total_SNPs"] += 1
        elif var_type == "MNP":
            result["Total_MNPs"] += 1
        elif var_type == "INDEL":
            result["Total_INDELs"] += 1

        allele_dict = {0: ref}
        for idx, alt in enumerate(alts, start=1):
            allele_dict[idx] = alt

        aligned_alleles = align_alleles_with_mafft(allele_dict)
        aln_len = len(next(iter(aligned_alleles.values())))

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
                if var_type == "SNP":
                    result[f"{s}_snps"] += 1
                elif var_type == "MNP":
                    result[f"{s}_mnps"] += 1
                elif var_type == "INDEL":
                    result[f"{s}_indels"] += 1

                if allele_counts.get(allele_idx, 0) == 1:
                    result[f"{s}_unique_alleles"] += 1

                seq_str = aligned_alleles.get(allele_idx, "-" * aln_len)
                sample_alignments[s].append(seq_str)
                result[f"{s}_bp"] += len(seq_str.replace("-", ""))

                var_gt_row.append(allele_idx)
                var_alt_row.append(1 if allele_idx > 0 else 0)
            else:
                seq_str = "-" * aln_len
                sample_alignments[s].append(seq_str)
                var_gt_row.append(np.nan)
                var_alt_row.append(np.nan)

        gt_matrix.append(var_gt_row)
        alt_bool_matrix.append(var_alt_row)

    if result["Total_variants"] == 0:
        return result

    # Total window sequence alignment
    full_sample_seqs = {s: "".join(sample_alignments[s]) for s in active_samples}
    result["Total_bp"] = len(next(iter(full_sample_seqs.values())))

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
        # 2 samples: construct simple unrooted Newick string (S1:d/2, S2:d/2) and assign patristic distance
        s1, s2 = active_samples[0], active_samples[1]
        d = seq_dist_matrix[0, 1]
        half_d = round(d / 2.0, 6)
        result["Tree_newick"] = f"({s1}:{half_d},{s2}:{half_d});"
        result[f"{s1}-{s2}_patristic"] = round(d, 6)

    return result
