import unittest
import os
import tempfile
import pysam
import pantreelib
from scripts import local


class TestPantreeLocal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.TemporaryDirectory()
        cls.vcf_path = os.path.join(cls.temp_dir.name, "test_local.vcf")
        cls.vcf_gz_path = os.path.join(cls.temp_dir.name, "test_local.vcf.gz")

        # Create sample VCF content:
        # pos 100: SAMPLE1=0/0, SAMPLE2=0/0, SAMPLE3=1/1 (SAMPLE1,SAMPLE2 share 0; SAMPLE3 shares 1)
        # pos 200: SAMPLE1=1/1, SAMPLE2=1/1, SAMPLE3=1/1 (All share 1 -> SAMPLE1,SAMPLE2,SAMPLE3)
        # pos 300: SAMPLE1=0/0, SAMPLE2=1/1, SAMPLE3=2/2 (Multiallelic: SAMPLE1 shares 0, SAMPLE2 shares 1, SAMPLE3 shares 2)
        # pos 400: SAMPLE1=./., SAMPLE2=0/0, SAMPLE3=0/1 (SAMPLE1 missing, SAMPLE2 shares 0, SAMPLE3 shares 1)
        # pos 500: SAMPLE1=0/0, SAMPLE2=0/0, SAMPLE3=0/0 (All share 0)
        vcf_content = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=1000>\n"
            "##contig=<ID=chr2,length=1000>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\tSAMPLE3\n"
            "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\n"
            "chr1\t200\t.\tC\tG\t.\tPASS\t.\tGT\t1/1\t1/1\t1/1\n"
            "chr1\t300\t.\tA\tC,G\t.\tPASS\t.\tGT\t0/0\t1/1\t2/2\n"
            "chr1\t400\t.\tC\tT\t.\tPASS\t.\tGT\t./.\t0/0\t0/1\n"
            "chr1\t500\t.\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\n"
        )
        with open(cls.vcf_path, "w") as f:
            f.write(vcf_content)

        pysam.tabix_compress(cls.vcf_path, cls.vcf_gz_path, force=True)
        pysam.tabix_index(cls.vcf_gz_path, preset="vcf", force=True)

    @classmethod
    def tearDownClass(cls):
        cls.temp_dir.cleanup()

    def test_process_local_window_basic(self):
        # Window chr1:1-300: 3 variants (100, 200, 300)
        # pos 100: {SAMPLE1, SAMPLE2}: 0, {SAMPLE3}: 1
        # pos 200: {SAMPLE1, SAMPLE2, SAMPLE3}: 1
        # pos 300: {SAMPLE1}: 0, {SAMPLE2}: 1, {SAMPLE3}: 2
        # Expected counts:
        # SAMPLE1,SAMPLE2 -> pos 100 (Count=1, Called=3)
        # SAMPLE3 -> pos 100, pos 300 (Count=2, Called=3)
        # SAMPLE1,SAMPLE2,SAMPLE3 -> pos 200 (Count=1, Called=3)
        # SAMPLE1 -> pos 300 (Count=1, Called=3)
        # SAMPLE2 -> pos 300 (Count=1, Called=3)
        task = (self.vcf_gz_path, "chr1", 1, 300, None, False)
        rows = pantreelib.process_local_window_task(task)

        row_map = {r["Group"]: r for r in rows}
        self.assertEqual(len(rows), 5)
        self.assertEqual(row_map["SAMPLE1,SAMPLE2"]["Count"], 1)
        self.assertEqual(row_map["SAMPLE1,SAMPLE2"]["Total_variants"], 3)
        self.assertEqual(row_map["SAMPLE1,SAMPLE2"]["Called_variants"], 3)

        self.assertEqual(row_map["SAMPLE3"]["Count"], 2)
        self.assertEqual(row_map["SAMPLE3"]["Called_variants"], 3)

        self.assertEqual(row_map["SAMPLE1,SAMPLE2,SAMPLE3"]["Count"], 1)
        self.assertEqual(row_map["SAMPLE1,SAMPLE2,SAMPLE3"]["Called_variants"], 3)

        self.assertEqual(row_map["SAMPLE1"]["Count"], 1)
        self.assertEqual(row_map["SAMPLE2"]["Count"], 1)

    def test_missing_data_and_called_variants(self):
        # Window chr1:350-450: pos 400 only
        # pos 400: SAMPLE1 missing (./.), SAMPLE2 has 0, SAMPLE3 has 1 (from 0/1)
        # Groups: SAMPLE2 (allele 0, count 1), SAMPLE3 (allele 1, count 1)
        # SAMPLE1 has no called variant at pos 400
        task = (self.vcf_gz_path, "chr1", 350, 450, None, False)
        rows = pantreelib.process_local_window_task(task)

        row_map = {r["Group"]: r for r in rows}
        self.assertEqual(len(rows), 2)
        self.assertIn("SAMPLE2", row_map)
        self.assertIn("SAMPLE3", row_map)
        self.assertEqual(row_map["SAMPLE2"]["Count"], 1)
        self.assertEqual(row_map["SAMPLE2"]["Total_variants"], 1)
        self.assertEqual(row_map["SAMPLE2"]["Called_variants"], 1)

        # Now include a wider window with pos 100, 200, 300, 400
        # Total_variants = 4
        # For SAMPLE1, called at 100, 200, 300 (not 400) -> Called_variants = 3
        # For group SAMPLE1,SAMPLE2 -> called at 100, 200, 300 (not 400) -> Called_variants = 3
        # For group SAMPLE2,SAMPLE3 -> called at 100, 200, 300, 400 -> Called_variants = 4
        task2 = (self.vcf_gz_path, "chr1", 1, 400, None, False)
        rows2 = pantreelib.process_local_window_task(task2)
        row_map2 = {r["Group"]: r for r in rows2}

        self.assertEqual(row_map2["SAMPLE1,SAMPLE2"]["Total_variants"], 4)
        self.assertEqual(row_map2["SAMPLE1,SAMPLE2"]["Called_variants"], 3)
        self.assertEqual(row_map2["SAMPLE2"]["Called_variants"], 4)

    def test_include_ref(self):
        # Window chr1:1-200 with include_ref=True
        # active_samples: REF, SAMPLE1, SAMPLE2, SAMPLE3
        # pos 100: REF=0, SAMPLE1=0, SAMPLE2=0, SAMPLE3=1
        # -> REF,SAMPLE1,SAMPLE2 share 0 (Count=1)
        # -> SAMPLE3 shares 1 (Count=1)
        # pos 200: REF=0, SAMPLE1=1, SAMPLE2=1, SAMPLE3=1
        # -> REF shares 0 (Count=1)
        # -> SAMPLE1,SAMPLE2,SAMPLE3 share 1 (Count=1)
        task = (self.vcf_gz_path, "chr1", 1, 200, None, True)
        rows = pantreelib.process_local_window_task(task)
        row_map = {r["Group"]: r for r in rows}

        self.assertIn("REF,SAMPLE1,SAMPLE2", row_map)
        self.assertIn("SAMPLE3", row_map)
        self.assertIn("REF", row_map)
        self.assertIn("SAMPLE1,SAMPLE2,SAMPLE3", row_map)
        self.assertEqual(row_map["REF,SAMPLE1,SAMPLE2"]["Count"], 1)
        self.assertEqual(row_map["SAMPLE1,SAMPLE2,SAMPLE3"]["Count"], 1)

    def test_samples_subset(self):
        # Subset SAMPLE1,SAMPLE3 only
        # pos 100: SAMPLE1=0, SAMPLE3=1 -> SAMPLE1 (0), SAMPLE3 (1)
        # pos 200: SAMPLE1=1, SAMPLE3=1 -> SAMPLE1,SAMPLE3 (1)
        task = (self.vcf_gz_path, "chr1", 1, 200, ["SAMPLE1", "SAMPLE3"], False)
        rows = pantreelib.process_local_window_task(task)
        row_map = {r["Group"]: r for r in rows}

        self.assertEqual(len(rows), 3)
        self.assertIn("SAMPLE1", row_map)
        self.assertIn("SAMPLE3", row_map)
        self.assertIn("SAMPLE1,SAMPLE3", row_map)

    def test_empty_window_omitted(self):
        # Window on chr2 has 0 variants
        task = (self.vcf_gz_path, "chr2", 1, 500, None, False)
        rows = pantreelib.process_local_window_task(task)
        self.assertEqual(rows, [])

    def test_local_run_integration(self):
        out_tsv = os.path.join(self.temp_dir.name, "out_local.tsv")

        class Args:
            input_vcf = self.vcf_gz_path
            output_tsv = out_tsv
            samples = None
            window_size = 300
            overlap = 0
            fixed_variants = False
            include_ref = False
            threads = 2
            verbose = False

        local.run(Args())
        self.assertTrue(os.path.exists(out_tsv))

        with open(out_tsv) as f:
            lines = [line.strip().split("\t") for line in f if line.strip()]

        header = lines[0]
        expected_header = [
            "Chromosome", "Window_start", "Window_end",
            "Total_variants", "Called_variants", "Group", "Count"
        ]
        self.assertEqual(header, expected_header)

        # Check rows exist and are sorted by Chromosome, Window_start, Group
        data_rows = lines[1:]
        self.assertGreater(len(data_rows), 0)
        sort_keys = [(r[0], int(r[1]), r[5]) for r in data_rows]
        self.assertEqual(sort_keys, sorted(sort_keys))

        # Check that empty chr2 windows were omitted
        chroms = {r[0] for r in data_rows}
        self.assertNotIn("chr2", chroms)


if __name__ == "__main__":
    unittest.main()
