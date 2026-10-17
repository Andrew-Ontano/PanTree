import unittest
import os
import tempfile
import pysam
import pantreelib
from scripts import calculate


class TestPantreeCalculate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.TemporaryDirectory()
        cls.vcf_path = os.path.join(cls.temp_dir.name, "test.vcf")
        cls.vcf_gz_path = os.path.join(cls.temp_dir.name, "test.vcf.gz")

        # Create sample VCF content
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

    @classmethod
    def tearDownClass(cls):
        cls.temp_dir.cleanup()

    def test_variant_classification(self):
        self.assertEqual(pantreelib.classify_variant("A", ["T"]), "SNP")
        self.assertEqual(pantreelib.classify_variant("AC", ["TG"]), "MNP")
        self.assertEqual(pantreelib.classify_variant("A", ["ATT", "G"]), "INDEL")

    def test_process_window_task(self):
        task = (self.vcf_gz_path, "chr1", 1, 500, None, False)
        result = pantreelib.process_window_task(task)

        self.assertEqual(result["Chromosome"], "chr1")
        self.assertEqual(result["Total_variants"], 4)
        self.assertEqual(result["Total_SNPs"], 2)
        self.assertEqual(result["Total_MNPs"], 1)
        self.assertEqual(result["Total_INDELs"], 1)
        self.assertEqual(result["SAMPLE1_called_sites"], 3)
        self.assertEqual(result["SAMPLE2_called_sites"], 4)
        self.assertEqual(result["SAMPLE3_called_sites"], 4)

    def test_calculate_run_integration(self):
        out_tsv = os.path.join(self.temp_dir.name, "out.tsv")

        class Args:
            input_vcf = self.vcf_gz_path
            output_tsv = out_tsv
            samples = "SAMPLE1,SAMPLE2"
            window_size = 500
            overlap = 0
            fixed_variants = False
            include_ref = True
            threads = 1
            verbose = False

        calculate.run(Args())
        self.assertTrue(os.path.exists(out_tsv))

        with open(out_tsv) as f:
            lines = [line.strip().split("\t") for line in f if line.strip()]

        header = lines[0]
        self.assertIn("Chromosome", header)
        self.assertIn("REF_called_sites", header)
        self.assertIn("SAMPLE1_called_sites", header)
        self.assertIn("SAMPLE2_called_sites", header)
        self.assertNotIn("SAMPLE3_called_sites", header)  # SAMPLE3 excluded via -s


if __name__ == "__main__":
    unittest.main()
