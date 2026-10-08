# by_label_null_partial (validator regression, not a generated fixture)

A real answer of `metagraph pattern --json` (the CLI, built at d0491f7d) on a tiny index the outside review GPT-2
built for its recheck: a graph `varying.dbg` with one 512 KiB label in a row-diff annotation (`long-label.row_diff
.annodbg`), `--pattern-build-mask --pattern-min-information-bits 0`, the request in `request.json`.

Its second pattern (`AAA`) shows what no fixture of the mini index can: in `partial`, the memory account could not
hold `by_label` (SPEC §14.4), so `by_label: null`, every context read `labels_status: "output_budget"`, `stop
{output, max_memory}`, and `counts.labels` still `exact`. `test_pattern_fixtures.py`
(`test_by_label_null_in_partial_is_valid`) checks that the validator accepts it. `scripts/traversal/pattern_fixtures.py`
neither writes nor checks this directory.
