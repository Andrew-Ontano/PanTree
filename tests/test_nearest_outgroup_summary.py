import unittest
import os
import tempfile
import pandas as pd
from utils.nearest_outgroup_summary import (
    calculate_nearest_outgroup_summary,
    plot_grouped_stacked_bar,
    get_pair_column
)


class TestNearestOutgroupSummary(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        # Synthetic DataFrame mimicking output of pantree.py calculate
        data = {
            "Chromosome": ["chr1", "chr1", "chr2"],
            "Window_start": [1, 501, 1],
            "Window_end": [500, 1000, 500],
            "Total_variants": [10, 20, 30],
            "SAMPLE1_called_sites": [10, 20, 30],
            "SAMPLE2_called_sites": [10, 20, 30],
            "OUT1_called_sites": [10, 20, 30],
            "OUT2_called_sites": [10, 20, 30],
            # Window 1 (chr1): SAMPLE1 closer to OUT1 (0.1 vs 0.5)
            # Window 2 (chr1): SAMPLE1 closer to OUT2 (0.4 vs 0.2)
            # Window 3 (chr2): SAMPLE1 tie with OUT1 & OUT2 (0.3 vs 0.3)
            "SAMPLE1-OUT1_manhattan": [0.1, 0.4, 0.3],
            "SAMPLE1-OUT2_manhattan": [0.5, 0.2, 0.3],
            "SAMPLE2-OUT1_manhattan": [0.2, 0.2, 0.4],
            "SAMPLE2-OUT2_manhattan": [0.3, 0.5, 0.1],
        }
        self.df = pd.DataFrame(data)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_get_pair_column(self):
        cols = self.df.columns
        self.assertEqual(get_pair_column("SAMPLE1", "OUT1", "manhattan", cols), "SAMPLE1-OUT1_manhattan")
        self.assertEqual(get_pair_column("OUT1", "SAMPLE1", "manhattan", cols), "SAMPLE1-OUT1_manhattan")
        self.assertIsNone(get_pair_column("SAMPLE1", "NONEXISTENT", "manhattan", cols))

    def test_calculate_nearest_outgroup_summary(self):
        selected = ["SAMPLE1"]
        non_selected = ["OUT1", "OUT2"]
        measures = ["manhattan"]

        summary = calculate_nearest_outgroup_summary(self.df, selected, non_selected, measures)

        # Scopes: chr1, chr2, PANGENOME
        self.assertIn("Scope", summary.columns)
        self.assertIn("Total_windows_closest", summary.columns)
        self.assertIn("Total_variants_closest", summary.columns)
        self.assertIn("Proportion_windows_closest", summary.columns)
        self.assertIn("Proportion_variants_closest", summary.columns)

        # chr1 checks for SAMPLE1:
        # Window 1 (10 vars): closest to OUT1 (win=1, vars=10)
        # Window 2 (20 vars): closest to OUT2 (win=1, vars=20)
        chr1_out1 = summary[(summary["Scope"] == "chr1") & (summary["Outgroup"] == "OUT1")].iloc[0]
        chr1_out2 = summary[(summary["Scope"] == "chr1") & (summary["Outgroup"] == "OUT2")].iloc[0]

        self.assertEqual(chr1_out1["Total_windows_closest"], 1.0)
        self.assertEqual(chr1_out1["Total_variants_closest"], 10.0)
        self.assertAlmostEqual(chr1_out1["Proportion_windows_closest"], 1.0 / 2.0, places=5)
        self.assertAlmostEqual(chr1_out1["Proportion_variants_closest"], 10.0 / 30.0, places=5)

        self.assertEqual(chr1_out2["Total_windows_closest"], 1.0)
        self.assertEqual(chr1_out2["Total_variants_closest"], 20.0)
        self.assertAlmostEqual(chr1_out2["Proportion_windows_closest"], 1.0 / 2.0, places=5)
        self.assertAlmostEqual(chr1_out2["Proportion_variants_closest"], 20.0 / 30.0, places=5)

        # chr2 checks for SAMPLE1:
        # Window 3 (30 vars): tie (0.3 vs 0.3) -> 0.5 win each, 15 vars each
        chr2_out1 = summary[(summary["Scope"] == "chr2") & (summary["Outgroup"] == "OUT1")].iloc[0]
        self.assertEqual(chr2_out1["Total_windows_closest"], 0.5)
        self.assertEqual(chr2_out1["Total_variants_closest"], 15.0)

        # PANGENOME checks:
        # Total windows: 3. Total variants: 60.
        # OUT1 total windows = 1.0 + 0.5 = 1.5; total variants = 10 + 15 = 25.
        pan_out1 = summary[(summary["Scope"] == "PANGENOME") & (summary["Outgroup"] == "OUT1")].iloc[0]
        self.assertEqual(pan_out1["Total_windows_closest"], 1.5)
        self.assertEqual(pan_out1["Total_variants_closest"], 25.0)
        self.assertAlmostEqual(pan_out1["Proportion_windows_closest"], 1.5 / 3.0, places=5)
        self.assertAlmostEqual(pan_out1["Proportion_variants_closest"], 25.0 / 60.0, places=5)

    def test_plot_generation(self):
        selected = ["SAMPLE1", "SAMPLE2"]
        non_selected = ["OUT1", "OUT2"]
        measures = ["manhattan"]

        summary = calculate_nearest_outgroup_summary(self.df, selected, non_selected, measures)
        plot_path = os.path.join(self.temp_dir.name, "test_plot.png")

        plot_grouped_stacked_bar(summary, "manhattan", selected, non_selected, plot_path)
        self.assertTrue(os.path.exists(plot_path))
        self.assertGreater(os.path.getsize(plot_path), 0)


if __name__ == "__main__":
    unittest.main()
