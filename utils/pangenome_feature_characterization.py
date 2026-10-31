import pysam
import gffutils
import argparse
import os
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description="Characterize pangenome VCF variants by GFF annotations.")
    parser.add_argument("-v", "--vcf", required=True, help="Input gzipped VCF (indexed)")
    parser.add_argument("-f", "--fasta", required=True, help="Reference FASTA (indexed)")
    parser.add_argument("-g", "--gff", required=True, help="GFF annotation file")
    parser.add_argument("-o", "--outdir", default="gene_alignments", help="Output directory for FASTA alignments")
    parser.add_argument("-a", "--align", action="store_true",
                        help="Preserve alignment gaps for indels to sync reading frames")
    parser.add_argument("-s", "--source", default="RefSeq",
                        help="GFF source filter to avoid duplicate models (default: RefSeq)")
    return parser.parse_args()


def init_gff_db(gff_path):
    db_path = gff_path + ".db"
    if not os.path.exists(db_path):
        gffutils.create_db(gff_path, db_path, force=True, keep_order=True, merge_strategy='merge',
                           sort_attribute_values=True)
    return gffutils.FeatureDB(db_path, keep_order=True)


def reconstruct_sequences(ref_seq, ref_start, variants, samples, align=False):
    """Reconstructs reference and sample sequences, optionally preserving gaps for alignment."""
    sequences = {s: [] for s in samples + ["Reference"]}
    curr_pos = ref_start

    for v in variants:
        v_start = v.pos - 1

        # Append conserved reference sequence between variants
        if v_start > curr_pos:
            chunk = ref_seq[curr_pos - ref_start: v_start - ref_start]
            for s in sequences:
                sequences[s].append(chunk)
        elif v_start < curr_pos:
            continue  # Skip overlapping variant representations

        # Determine max allele length for alignment padding
        max_len = len(v.ref)
        if align:
            for s in samples:
                gt = v.samples[s]['GT']
                if gt is not None and any(g > 0 for g in gt if g is not None):
                    alt_len = max([len(v.alts[g - 1]) for g in gt if g is not None and g > 0])
                    max_len = max(max_len, alt_len)

        # Process Reference Allele
        ref_allele = v.ref
        if align:
            ref_allele = ref_allele.ljust(max_len, '-')
        sequences["Reference"].append(ref_allele)

        # Process Sample Alleles
        for s in samples:
            gt = v.samples[s]['GT']
            if gt is None or not any(g > 0 for g in gt if g is not None):
                samp_allele = v.ref
            else:
                alt_idx = [g - 1 for g in gt if g is not None and g > 0][0]
                samp_allele = v.alts[alt_idx]

            if align:
                samp_allele = samp_allele.ljust(max_len, '-')
            sequences[s].append(samp_allele)

        curr_pos = v_start + len(v.ref)

    # Append remaining sequence after the last variant
    if curr_pos < ref_start + len(ref_seq):
        chunk = ref_seq[curr_pos - ref_start:]
        for s in sequences:
            sequences[s].append(chunk)

    return {s: "".join(sequences[s]) for s in sequences}


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    vcf_in = pysam.VariantFile(args.vcf)
    fasta_in = pysam.FastaFile(args.fasta)
    db = init_gff_db(args.gff)

    samples = list(vcf_in.header.samples)

    sample_occupancy = defaultdict(lambda: {"total_alleles": 0, "genes_with_change": 0})
    gene_variant_counts = defaultdict(lambda: defaultdict(int))

    genes = db.features_of_type('gene')

    for gene in genes:
        # Filter by annotation source to avoid doubling models
        if gene.source != args.source:
            continue

        chrom = gene.seqid
        start, end = gene.start, gene.end
        strand = gene.strand
        gene_id = gene.id or gene.attributes.get('Name', [f"{chrom}_{start}"])[0]

        try:
            records = list(vcf_in.fetch(chrom, start, end))
        except ValueError:
            continue

        if not records:
            continue

        ref_seq = fasta_in.fetch(chrom, start - 1, end)
        gene_has_change = defaultdict(bool)

        for rec in records:
            var_type = "SNP" if len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts) else "INDEL/SV"
            gene_variant_counts[gene_id][var_type] += 1

            for sample in samples:
                gt = rec.samples[sample]['GT']
                if gt is not None and any(g > 0 for g in gt if g is not None):
                    sample_occupancy[sample]["total_alleles"] += 1
                    gene_has_change[sample] = True

        for sample in samples:
            if gene_has_change[sample]:
                sample_occupancy[sample]["genes_with_change"] += 1

        # Reconstruct alignments with the joint function
        seq_dict = reconstruct_sequences(ref_seq, start - 1, records, samples, align=args.align)

        # Write gene alignment FASTA
        ref_description = f"{chrom};{start};{end};{strand}"
        alignment_records = [SeqRecord(Seq(seq_dict["Reference"]), id="Reference", description=ref_description)]

        for sample in samples:
            if gene_has_change[sample]:
                alignment_records.append(SeqRecord(Seq(seq_dict[sample]), id=sample, description=""))

        fasta_out = os.path.join(args.outdir, f"{gene_id}_alignment.fasta")
        SeqIO.write(alignment_records, fasta_out, "fasta")

    print("=== Sample Variant Occupancy ===")
    print("Sample\tTotal Alleles\tGenes with ≥1 Change")
    for sample in samples:
        print(f"{sample}\t{sample_occupancy[sample]['total_alleles']}\t{sample_occupancy[sample]['genes_with_change']}")

    print("\n=== Gene Variant Breakdown ===")
    print("GeneID\tSNP\tINDEL/SV")
    for gene_id, counts in gene_variant_counts.items():
        print(f"{gene_id}\t{counts['SNP']}\t{counts['INDEL/SV']}")


if __name__ == "__main__":
    main()