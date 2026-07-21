"""Read the sample manifest and generate MetaGraph differential-rule JSON."""

from __future__ import annotations

import csv
import json
import re
from collections import OrderedDict
from pathlib import Path
from typing import Any, Mapping


SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
COMPARISON_ID_RE = SAMPLE_ID_RE

ASSEMBLY_DEFAULTS: dict[str, Any] = {
    "count_kmers": True,
    "clean": False,
    "family_wise_error_rate": 0.05,
    "test_by_unitig": False,
    "test_type": "nbinom_exact",
    "assemble_shared": False,
    "min_count": 2,
    "min_recurrence": 1,
    "min_in_recurrence": 0,
    "min_out_recurrence": 0,
    "max_in_recurrence": 2**64 - 1,
    "max_out_recurrence": 2**64 - 1,
}

SUPPORTED_TESTS = {
    "notest",
    "poisson_binom",
    "poisson_exact",
    "nbinom_exact",
    "mwu",
    "cmh",
    "cmh_binary",
    "fisher_binary",
}


def _require_identifier(value: str, kind: str) -> None:
    if not SAMPLE_ID_RE.fullmatch(value):
        raise ValueError(
            f"Invalid {kind} {value!r}; use letters, digits, '.', '_' or '-'"
        )


def load_samples(path: str | Path) -> OrderedDict[str, dict[str, Any]]:
    """Load a TSV manifest, accumulating repeated read rows per sample."""
    path = Path(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sample", "reads"}.issubset(reader.fieldnames):
            raise ValueError(
                f"{path} must be tab-separated and contain 'sample' and 'reads' columns"
            )

        samples: OrderedDict[str, dict[str, Any]] = OrderedDict()
        metadata_columns = [c for c in reader.fieldnames if c not in {"sample", "reads"}]
        if not metadata_columns:
            raise ValueError(f"{path} must contain at least one metadata column")

        for line_number, row in enumerate(reader, start=2):
            sample = (row.get("sample") or "").strip()
            reads = (row.get("reads") or "").strip()
            if not sample and not reads and all(not (row.get(c) or "").strip() for c in metadata_columns):
                continue
            _require_identifier(sample, f"sample ID on line {line_number}")
            if not reads:
                raise ValueError(f"Missing reads path for sample {sample!r} on line {line_number}")

            metadata = {c: (row.get(c) or "").strip() for c in metadata_columns}
            if sample in samples:
                if samples[sample]["metadata"] != metadata:
                    raise ValueError(
                        f"Repeated sample {sample!r} has inconsistent metadata on line {line_number}"
                    )
                samples[sample]["reads"].append(reads)
            else:
                samples[sample] = {"reads": [reads], "metadata": metadata}

    if not samples:
        raise ValueError(f"No samples found in {path}")
    return samples


def merged_assembly_config(
    global_config: Mapping[str, Any], comparison_config: Mapping[str, Any]
) -> dict[str, Any]:
    merged = dict(ASSEMBLY_DEFAULTS)
    merged.update(global_config)
    merged.update(comparison_config.get("assembly", {}))
    validate_assembly_config(merged)
    return merged


def validate_assembly_config(config: Mapping[str, Any]) -> None:
    if config.get("count_kmers") is not True:
        raise ValueError("This workflow requires assembly.count_kmers: true")
    if config.get("test_type") not in SUPPORTED_TESTS:
        raise ValueError(
            f"Unsupported test_type {config.get('test_type')!r}; "
            f"choose one of {sorted(SUPPORTED_TESTS)}"
        )
    alpha = float(config.get("family_wise_error_rate", 0))
    if not 0 < alpha <= 1:
        raise ValueError("family_wise_error_rate must be in (0, 1]")
    for key in (
        "min_count",
        "min_recurrence",
        "min_in_recurrence",
        "min_out_recurrence",
        "max_in_recurrence",
        "max_out_recurrence",
    ):
        if int(config[key]) < 0:
            raise ValueError(f"{key} must be non-negative")
    if int(config["min_in_recurrence"]) > int(config["max_in_recurrence"]):
        raise ValueError("min_in_recurrence exceeds max_in_recurrence")
    if int(config["min_out_recurrence"]) > int(config["max_out_recurrence"]):
        raise ValueError("min_out_recurrence exceeds max_out_recurrence")


def comparison_samples(
    samples: Mapping[str, Mapping[str, Any]], comparison: Mapping[str, Any]
) -> tuple[list[str], list[str]]:
    column = str(comparison["metadata_column"])
    in_values = {str(v) for v in comparison["in_values"]}
    out_values = {str(v) for v in comparison["out_values"]}
    if in_values & out_values:
        raise ValueError(
            f"in_values and out_values overlap for metadata column {column!r}"
        )

    missing_column = [s for s, record in samples.items() if column not in record["metadata"]]
    if missing_column:
        raise ValueError(f"Metadata column {column!r} is absent from the sample manifest")

    in_samples = [
        sample for sample, record in samples.items() if record["metadata"][column] in in_values
    ]
    out_samples = [
        sample for sample, record in samples.items() if record["metadata"][column] in out_values
    ]
    if not in_samples or not out_samples:
        raise ValueError(
            f"Comparison on {column!r} selected {len(in_samples)} in and "
            f"{len(out_samples)} out samples"
        )
    if set(in_samples) & set(out_samples):
        raise ValueError("A sample cannot belong to both sides of a comparison")
    return in_samples, out_samples


def _comparison_component(value: str) -> str:
    component = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._-")
    return component or "value"


def generate_auto_comparisons(
    samples: Mapping[str, Mapping[str, Any]], settings: Mapping[str, Any]
) -> OrderedDict[str, dict[str, Any]]:
    """Generate safe one-vs-rest comparisons from categorical metadata."""
    generated: OrderedDict[str, dict[str, Any]] = OrderedDict()
    if not settings.get("enabled", False):
        return generated

    available_columns = list(next(iter(samples.values()))["metadata"])
    columns = list(settings.get("metadata_columns") or available_columns)
    unknown_columns = [column for column in columns if column not in available_columns]
    if unknown_columns:
        raise ValueError(
            "Unknown auto-comparison metadata columns: "
            + ", ".join(repr(column) for column in unknown_columns)
        )

    missing_values = {
        str(value).strip().casefold()
        for value in settings.get("missing_values", ["", "NA", "NaN", "null"])
    }
    minimum_samples = int(settings.get("minimum_samples_per_group", 1))
    if minimum_samples < 1:
        raise ValueError("minimum_samples_per_group must be at least one")
    assembly_overrides = dict(settings.get("assembly", {}))

    for column in columns:
        categories: OrderedDict[str, list[str]] = OrderedDict()
        for sample, record in samples.items():
            value = str(record["metadata"][column]).strip()
            if value.casefold() in missing_values:
                continue
            categories.setdefault(value, []).append(sample)

        for category, in_samples in categories.items():
            out_values = [value for value in categories if value != category]
            out_sample_count = sum(len(categories[value]) for value in out_values)
            if len(in_samples) < minimum_samples or out_sample_count < minimum_samples:
                continue

            comparison_name = "_".join(
                (
                    _comparison_component(column),
                    _comparison_component(category),
                    "vs_rest",
                )
            )
            if comparison_name in generated:
                raise ValueError(
                    f"Auto-comparison name collision for {comparison_name!r}; "
                    "rename the metadata column or category"
                )

            comparison: dict[str, Any] = {
                "metadata_column": column,
                "in_values": [category],
                "out_values": out_values,
            }
            if assembly_overrides:
                comparison["assembly"] = dict(assembly_overrides)
            generated[comparison_name] = comparison

    return generated


def make_diff_document(
    comparison_name: str,
    in_labels: list[str],
    out_labels: list[str],
    assembly_config: Mapping[str, Any],
) -> dict[str, Any]:
    _require_identifier(comparison_name, "comparison name")
    validate_assembly_config(assembly_config)
    experiment = {key: assembly_config[key] for key in ASSEMBLY_DEFAULTS}
    experiment.update(
        {
            "name": comparison_name,
            "in": in_labels,
            "out": out_labels,
        }
    )
    return {"groups": [{"experiments": [experiment]}]}


def write_diff_document(document: Mapping[str, Any], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(document, handle, indent=2, sort_keys=True)
        handle.write("\n")
