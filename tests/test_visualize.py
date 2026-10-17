import unittest
import os
import tempfile
import pysam
import pandas as pd
from scripts import calculate, visualize
import pantreelib


class TestPantreeVisualize(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.TemporaryDirectory()
        cls.vcf_path = os.path.join(cls.temp_dir.name, "test.vcf")
        cls.vcf_gz_path = os.path.join(cls.temp_dir.name, "test.vcf.gz")
        cls.tsv_path = os.path.join(cls.temp_dir.name, "test_out.tsv")
        cls.plots_dir = os.path.join(cls.temp_dir.name, "plots")

        vcf_content = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=1000>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\tSAMPLE3\n"
            "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n"
            "chr1\t200\t.\tAC\tTG\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\n"
            "chr1\t300\t.\tA\tATT,G\t.\tPASS\t.\tGT\t0/0\t1/1\t2/2\n"
            "chr1\t400\t.\tC\tT\t.\tPASS\t.\tGT\t./.\t0/0\t0/1\n"
        )
        with open(cls.vcf_path, "w") as f:
            f.write(vcf_content)

        pysam.tabix_compress(cls.vcf_path, cls.vcf_gz_path, force=True)
        pysam.tabix_index(cls.vcf_gz_path, preset="vcf", force=True)

        class CalcArgs:
            input_vcf = cls.vcf_gz_path
            output_tsv = cls.tsv_path
            samples = None
            window_size = 250
            overlap = 0
            fixed_variants = False
            include_ref = True
            threads = 1
            verbose = False

        calculate.run(CalcArgs())

    @classmethod
    def tearDownClass(cls):
        cls.temp_dir.cleanup()

    def test_jaccard_non_ref_pairs_when_ref_included(self):
        df = pd.read_csv(self.tsv_path, sep='\t')
        row = df[df['Total_variants'] > 0].iloc[0]

        # REF-SAMPLE Jaccard should be NA (parsed as NaN by pandas)
        self.assertTrue(pd.isna(row['REF-SAMPLE1_jaccard']))

        # Non-REF sample pairs should have valid numerical Jaccard values
        self.assertFalse(pd.isna(row['SAMPLE1-SAMPLE2_jaccard']))
        self.assertFalse(pd.isna(row['SAMPLE1-SAMPLE3_jaccard']))

    def test_visualize_run(self):
        class VisArgs:
            input_tsv = self.tsv_path
            outdir = self.plots_dir
            samples = None
            focus = None
            frac = 0.1
            format = "png"
            verbose = False

        visualize.run(VisArgs())
        self.assertTrue(os.path.exists(self.plots_dir))

        files = os.listdir(self.plots_dir)
        expected_files = [
            "chr1_patristic_density_variants.png",
            "chr1_patristic_density_bp.png",
            "chr1_jaccard_density_variants.png",
            "chr1_jaccard_density_bp.png",
            "chr1_manhattan_density_variants.png",
            "chr1_manhattan_density_bp.png",
            "chr1_variant_occupancy.png",
            "chr1_unique_variants_proportion.png"
        ]

        for expected in expected_files:
            self.assertIn(expected, files)


if __name__ == "__main__":
    unittest.main()
