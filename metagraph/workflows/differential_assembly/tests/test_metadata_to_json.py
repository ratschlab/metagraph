import json
import tempfile
import unittest
from collections import OrderedDict
from pathlib import Path

from workflow.scripts.metadata_to_json import (
    comparison_samples,
    generate_auto_comparisons,
    load_samples,
    make_diff_document,
    merged_assembly_config,
    write_diff_document,
)


class MetadataToJsonTest(unittest.TestCase):
    def test_manifest_comparison_and_document(self):
        samples = load_samples(Path(__file__).with_name("samples.tsv"))
        comparison = {
            "metadata_column": "phenotype",
            "in_values": ["case"],
            "out_values": ["control"],
        }
        in_samples, out_samples = comparison_samples(samples, comparison)
        self.assertEqual(in_samples, ["case_1", "case_2"])
        self.assertEqual(out_samples, ["control_1", "control_2"])

        assembly = merged_assembly_config(
            {"count_kmers": True, "test_type": "poisson_exact"}, comparison
        )
        document = make_diff_document(
            "case_vs_control",
            [f"results/{sample}.fasta.gz" for sample in in_samples],
            [f"results/{sample}.fasta.gz" for sample in out_samples],
            assembly,
        )
        experiment = document["groups"][0]["experiments"][0]
        self.assertTrue(experiment["count_kmers"])
        self.assertEqual(experiment["test_type"], "poisson_exact")
        self.assertEqual(len(experiment["in"]), 2)
        self.assertEqual(len(experiment["out"]), 2)

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "rules.json"
            write_diff_document(document, output)
            self.assertEqual(json.loads(output.read_text()), document)

    def test_auto_comparisons_are_safe_one_vs_rest(self):
        samples = OrderedDict(
            (
                ("case_1", {"metadata": {"phenotype": "case group"}}),
                ("case_2", {"metadata": {"phenotype": "case group"}}),
                ("control_1", {"metadata": {"phenotype": "control/A"}}),
                ("control_2", {"metadata": {"phenotype": "control/A"}}),
                ("other_1", {"metadata": {"phenotype": "other"}}),
                ("missing_1", {"metadata": {"phenotype": "NA"}}),
            )
        )
        comparisons = generate_auto_comparisons(
            samples,
            {
                "enabled": True,
                "metadata_columns": ["phenotype"],
                "minimum_samples_per_group": 2,
                "assembly": {"test_by_unitig": True},
            },
        )

        self.assertEqual(
            list(comparisons),
            ["phenotype_case_group_vs_rest", "phenotype_control_A_vs_rest"],
        )
        case = comparisons["phenotype_case_group_vs_rest"]
        self.assertEqual(case["in_values"], ["case group"])
        self.assertEqual(case["out_values"], ["control/A", "other"])
        self.assertTrue(case["assembly"]["test_by_unitig"])


if __name__ == "__main__":
    unittest.main()
