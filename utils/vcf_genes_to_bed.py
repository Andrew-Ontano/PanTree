import pysam
import gffutils
import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Extract gene-localized pangenome variants to BED for JBrowse 2.")
    parser.add_argument("-v", "--vcf", required=True, help="Input gzipped VCF (indexed)")
    parser.add_argument("-g", "--gff", required=True, help="GFF annotation file")
    parser.add_argument("-o", "--outbed", default="gene_variants.bed", help="Output BED file")
    parser.add_argument("-s", "--source", default="RefSeq", help="GFF source filter (default: RefSeq)")
    return parser.parse_args()


def init_gff_db(gff_path):
    db_path = gff_path + ".db"
    if not os.path.exists(db_path):
        gffutils.create_db(gff_path, db_path, force=True, keep_order=True, merge_strategy='merge',
                           sort_attribute_values=True)
    return gffutils.FeatureDB(db_path, keep_order=True)


def main():
    args = parse_args()

    vcf_in = pysam.VariantFile(args.vcf)
    db = init_gff_db(args.gff)

    genes = db.features_of_type('gene')

    with open(args.outbed, 'w') as bed_out:
        for gene in genes:
            if gene.source != args.source:
                continue

            chrom = gene.seqid
            start, end = gene.start, gene.end
            strand = gene.strand
            gene_id = gene.id or gene.attributes.get('Name', [f"{chrom}_{start}"])[0]

            try:
                records = list(vcf_in.fetch(chrom, start, end))
            except ValueError:
                continue  # Chromosome not in VCF

            for rec in records:
                # Determine variant type and length for BED coordinates
                is_snp = len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts)
                var_type = "SNP" if is_snp else "INDEL/SV"

                # VCF is 1-based, BED is 0-based
                v_start = rec.pos - 1
                v_end = v_start + len(rec.ref)

                # Format: chrom, start, end, name, score (0), strand
                feature_name = f"{gene_id}_{var_type}"
                bed_out.write(f"{chrom}\t{v_start}\t{v_end}\t{feature_name}\t0\t{strand}\n")


if __name__ == "__main__":
    main()