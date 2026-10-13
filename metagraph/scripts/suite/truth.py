"""Ground truth of the suite: the genomes, what the index holds of each record, and the record sequences fetched
at their indexed accession.version from NCBI E-utilities (cached on disk with their sha256).

Inputs (all in scripts/suite/truth/, committed):
- search_service/genomes.json, manifest_v1.json: the search service's validation genomes and verified proteins
  (pinned copies, source commit in SOURCES below);
- index/genome_records.tsv, genome_columns.tsv, last_record_scan.jsonl: the records' column, seq_id and column
  sizes in refseq33m's `.seqs`, and the 2026-10-09 scan of last records (from the index defect report);
- cache/ (not committed): the fetched FASTA, one file per record (or region).

A record over WHOLE_RECORD_MAX_BP is fetched as the region genomes.json names (1-based inclusive), so that the
cache stays small; every coordinate the suite states is a coordinate of the whole record.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Dict, List, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH_DIR = os.path.join(HERE, "truth")
CACHE_DIR = os.environ.get("SUITE_TRUTH_CACHE", os.path.join(TRUTH_DIR, "cache"))

SOURCES = {
    "genomes": "metagraph-search-service scripts/pattern-validation/genomes.json @ "
               "20343460e658e1f88478227f79f19e1862b27cae",
    "manifest_v1": "metagraph-search-service scripts/pattern-validation/manifest/manifest_v1.json @ "
                   "20343460e658e1f88478227f79f19e1862b27cae (content_sha256 cfa1e5bd261e3531...)",
    "index_records": "refseq33m annotation.relaxed.relabeled.seqs, decoded 2026-10-09 "
                     "(~/.cache/metagraph-index-defects/refseq33m-NC_000932/data/records.tsv, last_records.tsv)",
    "last_record_scan": "the same report's scan_last_NC_multi.jsonl (first 200 bp of each last record vs staging)",
}

WHOLE_RECORD_MAX_BP = 10_000_000
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
K = 31


@dataclass
class Record:
    key: str
    label: str
    accession: str
    kingdom: str
    compartment: str
    gc_class: str
    transl_table: int
    column: str                      # the index column label (a taxid string)
    column_index: int
    seq_id: int                      # the record's index in its column (.seqs)
    records_in_column: int
    length: int                      # bp (the .seqs k-mers + k - 1)
    last_in_column: bool
    last_record_scan: Optional[str]  # status of the 2026-10-09 scan for a last record ("ok", "header_none")
    region: Optional[List[int]]      # [start, end] 1-based inclusive when only a region is fetched
    genomes_json_taxid: Optional[int] = None
    seq: str = ""                    # the fetched sequence (whole record or region), upper case
    offset: int = 0                  # record coordinate (0-based) of seq[0]
    sha256: str = ""
    fetch_error: str = ""

    @property
    def usable_for_positives(self) -> bool:
        """A record whose k-mers are annotated under its column: not a last record the scan found missing."""
        return not self.fetch_error and self.last_record_scan != "header_none"

    def describe(self) -> dict:
        return {
            "key": self.key, "accession": self.accession, "label": self.label, "taxid": self.column,
            "column_index": self.column_index, "seq_id": self.seq_id,
            "records_in_column": self.records_in_column, "length": self.length,
            "last_in_column": self.last_in_column, "last_record_scan": self.last_record_scan,
            "transl_table": self.transl_table, "kingdom": self.kingdom, "compartment": self.compartment,
            "gc_class": self.gc_class,
            "fetched": ({"region": self.region} if self.region else {"region": None}),
            "fetched_bp": len(self.seq), "sha256": self.sha256, "fetched_from": "efetch",
            **({"fetch_error": self.fetch_error} if self.fetch_error else {}),
            **({"genomes_json_taxid": self.genomes_json_taxid}
               if self.genomes_json_taxid is not None and str(self.genomes_json_taxid) != self.column else {}),
        }


def _read_tsv(path: str) -> List[dict]:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_genomes() -> List[dict]:
    with open(os.path.join(TRUTH_DIR, "search_service", "genomes.json")) as f:
        return json.load(f)["genomes"]


def load_proteins() -> List[dict]:
    with open(os.path.join(TRUTH_DIR, "search_service", "manifest_v1.json")) as f:
        return json.load(f)["proteins"]


def load_records() -> List[Record]:
    """The genomes with the index facts of their records (no sequence)."""
    recs = {r["header"]: r for r in _read_tsv(os.path.join(TRUTH_DIR, "index", "genome_records.tsv"))}
    cols = {r["column_index"]: r for r in _read_tsv(os.path.join(TRUTH_DIR, "index", "genome_columns.tsv"))}
    scan = {}
    with open(os.path.join(TRUTH_DIR, "index", "last_record_scan.jsonl")) as f:
        for line in f:
            j = json.loads(line)
            scan[j["header"]] = j["status"]
    out = []
    for g in load_genomes():
        r = recs[g["accession"]]
        c = cols[r["column_index"]]
        length = int(r["kmers"]) + K - 1
        last = c["last_header"] == g["accession"]
        n = int(c["records_in_column"])
        region = None
        if length > WHOLE_RECORD_MAX_BP:
            region = list(g["region"])
        out.append(Record(
            key=g["key"], label=g["label"], accession=g["accession"], kingdom=g["kingdom"],
            compartment=g["compartment"], gc_class=g["gc_class"], transl_table=int(g["transl_table"]),
            column=r["column_label"], column_index=int(r["column_index"]), seq_id=int(r["seq_id"]),
            records_in_column=n, length=length, last_in_column=last,
            last_record_scan=(scan.get(g["accession"]) if last and n > 1 else None),
            region=region, genomes_json_taxid=g.get("taxid")))
    return out


def cache_path(rec: Record) -> str:
    name = rec.accession if not rec.region else "%s_%d-%d" % (rec.accession, rec.region[0], rec.region[1])
    return os.path.join(CACHE_DIR, name + ".fa")


def _parse_fasta(text: str) -> (str, str):
    lines = text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("not a FASTA answer: %r" % text[:200])
    if sum(1 for l in lines if l.startswith(">")) != 1:
        raise ValueError("more than one FASTA record")
    return lines[0][1:], "".join(l.strip() for l in lines[1:]).upper()


def fetch_record(rec: Record, pause_s: float = 0.4, retries: int = 3) -> None:
    """efetch the record (or its region) unless cached; fills seq, offset, sha256 (or fetch_error)."""
    path = cache_path(rec)
    if not os.path.exists(path):
        params = {"db": "nuccore", "id": rec.accession, "rettype": "fasta", "retmode": "text",
                  "tool": "metagraph-suite"}
        if rec.region:
            params["seq_start"] = str(rec.region[0])
            params["seq_stop"] = str(rec.region[1])
        url = EUTILS + "?" + urllib.parse.urlencode(params)
        last_err = ""
        for attempt in range(retries):
            time.sleep(pause_s * (1 + 2 * attempt))   # at most ~2.5 requests per second
            try:
                with urllib.request.urlopen(url, timeout=300) as r:
                    text = r.read().decode("ascii")
                header, seq = _parse_fasta(text)
                if not header.startswith(rec.accession):
                    raise ValueError("answer is %r, not %s" % (header[:60], rec.accession))
                os.makedirs(CACHE_DIR, exist_ok=True)
                tmp = path + ".tmp"
                with open(tmp, "w") as f:
                    f.write(">" + header + "\n" + seq + "\n")
                os.replace(tmp, path)
                last_err = ""
                break
            except (urllib.error.URLError, ValueError, OSError) as e:
                last_err = "%s: %s" % (type(e).__name__, e)
        if last_err:
            rec.fetch_error = last_err
            return
    with open(path) as f:
        _, seq = _parse_fasta(f.read())
    expected = (rec.region[1] - rec.region[0] + 1) if rec.region else rec.length
    if len(seq) != expected:
        rec.fetch_error = "fetched %d bp, the index states %d" % (len(seq), expected)
        return
    rec.seq = seq
    rec.offset = (rec.region[0] - 1) if rec.region else 0
    rec.sha256 = hashlib.sha256(seq.encode()).hexdigest()


def load_with_sequences(fetch: bool = True, keys: Optional[List[str]] = None) -> List[Record]:
    recs = load_records()
    if keys:
        recs = [r for r in recs if r.key in keys]
    for r in recs:
        if fetch or os.path.exists(cache_path(r)):
            fetch_record(r)
        else:
            r.fetch_error = "not cached (run with --fetch)"
    return recs


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="fetch and cache the suite's record sequences")
    ap.add_argument("--keys", help="comma-separated genome keys (default: all)")
    args = ap.parse_args(argv)
    keys = args.keys.split(",") if args.keys else None
    recs = load_with_sequences(fetch=True, keys=keys)
    for r in recs:
        print("%-5s %-15s %10d bp %s %s" % (r.key, r.accession, len(r.seq),
                                            (r.sha256[:12] or "-"), r.fetch_error), flush=True)
    return 1 if any(r.fetch_error for r in recs) else 0


if __name__ == "__main__":
    sys.exit(main())
