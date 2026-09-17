import unittest
import os
import tempfile
import pandas as pd
import pysam
import sys

# Add repo root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from PanTree.utils.vcf_dist import main as vcf_dist_main


def create_test_vcf(vcf_path, num_variants=25, chrom_length=100000):
    vcf_text = (
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID=chr1,length={chrom_length}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3\n"
    )
    # Generate variant positions spaced out
    positions = [100 + i * 2000 for i in range(num_variants)]
    for pos in positions:
        vcf_text += f"chr1\t{pos}\t.\tA\tG\t100\tPASS\t.\tGT\t0/0\t0/1\t1/1\n"

    uncompressed = vcf_path.replace(".vcf.gz", ".vcf")
    with open(uncompressed, "w") as f:
        f.write(vcf_text)

    pysam.tabix_compress(uncompressed, vcf_path, force=True)
    pysam.tabix_index(vcf_path, preset="vcf", force=True)
    if os.path.exists(uncompressed):
        os.remove(uncompressed)


class TestVcfDist(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.vcf_path = os.path.join(self.temp_dir.name, "test.vcf.gz")
        create_test_vcf(self.vcf_path, num_variants=25, chrom_length=50000)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_genomic_window_no_overlap(self):
        output_tsv = os.path.join(self.temp_dir.name, "out1.tsv")
        sys.argv = [
            "vcf_dist.py",
            "-v", self.vcf_path,
            "-w", "10000",
            "-o", output_tsv,
            "-t", "1"
        ]
        vcf_dist_main()
        df = pd.read_csv(output_tsv, sep="\t")
        # 50000 / 10000 = 5 windows
        self.assertEqual(len(df), 5)
        self.assertEqual(df.iloc[0]["Window_start"], 1)
        self.assertEqual(df.iloc[0]["Window_end"], 10000)
        self.assertEqual(df.iloc[1]["Window_start"], 10001)

    def test_genomic_window_with_overlap(self):
        output_tsv = os.path.join(self.temp_dir.name, "out2.tsv")
        sys.argv = [
            "vcf_dist.py",
            "-v", self.vcf_path,
            "-w", "10000",
            "-l", "2000",
            "-o", output_tsv,
            "-t", "1"
        ]
        vcf_dist_main()
        df = pd.read_csv(output_tsv, sep="\t")
        # Step is 2000: start = 1, 2001, 4001, ..., up to <= 50000
        starts = df["Window_start"].tolist()
        self.assertEqual(starts[0], 1)
        self.assertEqual(starts[1], 2001)
        self.assertEqual(starts[2], 4001)
        self.assertEqual(df.iloc[0]["Window_end"], 10000)
        self.assertEqual(df.iloc[1]["Window_end"], 12000)

    def test_fixed_variant_window_no_overlap(self):
        output_tsv = os.path.join(self.temp_dir.name, "out3.tsv")
        sys.argv = [
            "vcf_dist.py",
            "-v", self.vcf_path,
            "-w", "10",
            "-f",
            "-o", output_tsv,
            "-t", "1"
        ]
        vcf_dist_main()
        df = pd.read_csv(output_tsv, sep="\t")
        # 25 variants with window size 10 -> windows for variants 0-9, 10-19, 20-24
        self.assertEqual(len(df), 3)
        # Variant 0 position: 100, variant 9 position: 100 + 9*2000 = 18100
        self.assertEqual(df.iloc[0]["Window_start"], 100)
        self.assertEqual(df.iloc[0]["Window_end"], 18100)
        self.assertEqual(df.iloc[0]["Variant_count"], 10)

        # Second window: variant 10 position: 100 + 10*2000 = 20100
        self.assertEqual(df.iloc[1]["Window_start"], 20100)
        self.assertEqual(df.iloc[1]["Variant_count"], 10)

        # Third partial window: 5 variants
        self.assertEqual(df.iloc[2]["Variant_count"], 5)

    def test_fixed_variant_window_with_overlap(self):
        output_tsv = os.path.join(self.temp_dir.name, "out4.tsv")
        sys.argv = [
            "vcf_dist.py",
            "-v", self.vcf_path,
            "-w", "10",
            "-l", "2",
            "-f",
            "-o", output_tsv,
            "-t", "1"
        ]
        vcf_dist_main()
        df = pd.read_csv(output_tsv, sep="\t")
        # Step = 2, variants = 25
        # i = 0 (var 0..9), i = 2 (var 2..11), i = 4 (var 4..13), etc.
        self.assertEqual(df.iloc[0]["Window_start"], 100)  # var 0
        self.assertEqual(df.iloc[1]["Window_start"], 100 + 2 * 2000)  # var 2


if __name__ == "__main__":
    unittest.main()
