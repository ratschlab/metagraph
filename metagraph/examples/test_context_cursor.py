"""Exercise the JSONL example using freshly built coordinate indexes."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    cursor, metagraph = map(lambda p: str(Path(p).resolve()), sys.argv[1:])
    with tempfile.TemporaryDirectory(prefix="context-example-") as directory:
        root = Path(directory)
        fasta = root / "records.fa"
        fasta.write_text(">left_C\nCACGAG\n>left_T\nTACGAA\n>gapped\nAACGNTTTA\n>repeat\nAAAAAAA\n")
        graph, annotation = str(root / "graph"), str(root / "annotation")

        def build(*args):
            result = subprocess.run([metagraph, *map(str, args)], capture_output=True, text=True)
            assert result.returncode == 0, result.stderr

        build("build", "--in-ram", "--mask-dummy", "--mode", "basic", "--state", "stat",
              "-p", "1", "-k", "3", "-o", graph, fasta)
        build("annotate", "--anno-filename", "--coordinates", "-p", "1", "-i", graph + ".dbg",
              "-o", annotation, fasta)
        build("transform_anno", "--anno-type", "column_coord", "-p", "1", "-o", annotation,
              annotation + ".column.annodbg")
        build("annotate", "--anno-filename", "--index-header-coords", "-p", "1",
              "-i", graph + ".dbg", "-o", annotation, fasta)
        build("transform_anno", "--anno-type", "brwt_coord", "-p", "1", "-o", annotation,
              annotation + ".column.annodbg")

        valid = [
            {"seed": "ACGA", "steps": [{"side": "left", "base": "C"}]},
            {"seed": "ACGA", "radius": 0, "steps": [
                {"side": "left", "base": "c"}, {"side": "right", "base": "G", "radius": 10}]},
            {"seed": "ACGA", "steps": [{"side": "right", "base": "C"}]},
            {"seed": "ACG", "steps": [{"side": "left", "base": "A"}]},
            {"seed": "GGG"},
            {"seed": "AAA", "max_occurrences": 5},
        ]
        invalid = [
            {}, {"seed": 123}, {"seed": "AC"}, {"seed": "ACN"},
            {"seed": "A" * 10001}, {"seed": "AAA", "max_occurrences": 4},
            {"seed": "ACGA", "max_occurrences": 0}, {"seed": "ACGA", "max_occurrences": 100001},
            {"seed": "ACGA", "radius": -1}, {"seed": "ACGA", "radius": 10001},
            {"seed": "ACGA", "radius": 1.5}, {"seed": "ACGA", "radius": True},
            {"seed": "ACGA", "raduis": 10}, {"seed": "ACGA", "columns": "file"},
            {"seed": "ACGA", "columns": [1]}, {"seed": "ACGA", "columns": ["absent"]},
            {"seed": "ACGA", "steps": {}}, {"seed": "ACGA", "steps": [{}] * 101},
            {"seed": "ACGA", "steps": [{"side": "up", "base": "C"}]},
            {"seed": "ACGA", "steps": [{"side": "left", "base": "N"}]},
            {"seed": "ACGA", "steps": [{"side": "left", "base": "CC"}]},
            {"seed": "ACGA", "steps": [{"side": "right", "base": "C", "radius": -1}]},
            {"seed": "ACGA", "steps": [{"side": "left", "base": "C", "raduis": 1}]},
            {"seed": "ACGA", "steps": [False]},
        ]
        malformed = ['{"seed":"ACGA","seed":"AAA"}', '{"seed":"ACGA"} trailing',
                     '{"seed":"ACGA",}', '{"seed":', '[]']
        requests = root / "requests.jsonl"

        def run(kind, mmap, lines, expected_code):
            requests.write_text("\n".join(lines) + "\n")
            result = subprocess.run([cursor, "query", *( ["--mmap"] if mmap else []),
                "-i", graph + ".dbg", "-a", annotation + "." + kind + ".annodbg", str(requests)],
                capture_output=True, text=True)
            assert result.returncode == expected_code, result.stderr + result.stdout
            values = [json.loads(line) for line in result.stdout.splitlines()]
            assert len(values) == len(lines), values
            return values

        for kind in ("column_coord", "brwt_coord"):
            for mmap in (False, True):
                results = run(kind, mmap, list(map(json.dumps, valid)), 0)
                first = results[0]
                assert first["sequence"] == "ACGA" and len(first["occurrences"]) == 2
                assert first["closure"] == dict(left_added=0, right_added=0, left_stop="branch", right_stop="branch")
                chosen = first["steps"][0]
                assert chosen["sequence"] == "CACGAG" and chosen["seed_offset"] == 1
                assert chosen["occurrences"][0]["header"] == "left_C"
                assert chosen["closure"] == dict(left_added=0, right_added=1, left_stop="record_end", right_stop="record_end")
                limited = results[1]
                assert limited["closure"]["left_stop"] == "radius_limit"
                assert limited["steps"][0]["sequence"] == "CACGA"  # inherits radius=0
                assert limited["steps"][0]["closure"]["right_stop"] == "radius_limit"
                assert limited["steps"][1]["sequence"] == "CACGAG"
                failed = results[2]["steps"][0]
                assert not failed["extended"] and failed["closure"] is None
                assert failed["sequence"] == "ACGA" and failed["occurrences"] == results[2]["occurrences"]
                gapped = results[3]["steps"][0]
                assert gapped["sequence"] == "AACG"
                assert gapped["closure"]["right_stop"] == "unindexed"
                assert results[4]["occurrences"] == []
                assert results[4]["closure"]["left_stop"] == "no_occurrences"
                assert len(results[5]["occurrences"]) == 5
                assert results[5]["closure"]["left_stop"] == "record_end"

                lines = list(map(json.dumps, invalid)) + malformed + [json.dumps(valid[0])]
                errors = run(kind, mmap, lines, 1)
                assert all(isinstance(value.get("error"), str) for value in errors[:-1]), errors
                assert errors[-1] == first  # errors do not prevent later requests
        print(f"Context example: 4 backend/load combinations, {4 * (len(valid) + len(invalid) + len(malformed) + 1)} requests passed")


if __name__ == "__main__":
    main()
