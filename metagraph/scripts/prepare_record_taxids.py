#!/usr/bin/env python3
"""Join exact indexed headers to an NCBI RefSeq release catalogue, with complete coverage.

Export each chunk's inventory first with:
  metagraph stats --print-sequence-headers annotation.seqs > headers.tsv

This program creates a new output directory and uses a disk-backed SQLite join. It
never downloads data, edits an index, overwrites an output, or publishes a partial
mapping. Keep the catalogue from the SAME release as the indexed sequences.
"""

import argparse
import gzip
import hashlib
import json
import re
import sqlite3
from pathlib import Path


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def number(text, positive=False):
    if not re.fullmatch(r'0|[1-9][0-9]*', text):
        raise ValueError('Expected a canonical nonnegative integer')
    value = int(text)
    if (positive and not value) or value > 0xFFFFFFFF:
        raise ValueError('Integer outside the supported 32-bit range')
    return value


def prepare(catalogue, chunks, output):
    catalogue, output = Path(catalogue), Path(output)
    if not chunks:
        raise ValueError('At least one chunk is required')
    if any(not re.fullmatch(r'[A-Za-z0-9_-]+', name) for name in chunks):
        raise ValueError('Chunk names must contain only letters, digits, underscores or hyphens')
    output.mkdir(parents=False, exist_ok=False)
    report = {'format': 1, 'chunks_total': len(chunks), 'catalogue_sha256': sha256(catalogue),
              'catalogue_file': catalogue.name, 'chunks': {}, 'complete': False}
    print(json.dumps({'event': 'begin', 'chunks_total': len(chunks)}), flush=True)
    with sqlite3.connect(output / 'join.sqlite3') as db:
        db.executescript('''
            PRAGMA cache_size=-262144;
            PRAGMA temp_store=FILE;
            CREATE TABLE accession (name TEXT PRIMARY KEY, taxid INTEGER) WITHOUT ROWID;
            CREATE TABLE record (chunk TEXT, col INTEGER, seq INTEGER, accession TEXT,
                                 PRIMARY KEY (chunk, col, seq)) WITHOUT ROWID;
            CREATE TRIGGER consistent_taxid BEFORE UPDATE OF taxid ON accession
              WHEN OLD.taxid IS NOT NULL AND OLD.taxid != NEW.taxid
              BEGIN SELECT RAISE(ABORT, 'Conflicting taxids for an indexed accession'); END;
        ''')
        prefixes = set()
        for name, inventory in chunks.items():
            inventory = Path(inventory)
            count, previous = 0, None
            with inventory.open(encoding='utf-8', newline='') as stream:
                if stream.readline() != 'METAGRAPH_SEQUENCE_HEADERS_V1\n':
                    raise ValueError(f'{name}: unsupported header inventory')
                batch = []
                for line in stream:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) != 3 or not fields[2] or '\r' in line:
                        raise ValueError(f'{name}: invalid header inventory row')
                    col, seq = number(fields[0]), number(fields[1])
                    if previous is None:
                        valid = seq == 0
                    else:
                        valid = (col == previous[0] and seq == previous[1] + 1
                                 or col > previous[0] and seq == 0)
                    if not valid:
                        raise ValueError(f'{name}: unordered or missing record ids')
                    previous = (col, seq)
                    accession = fields[2]
                    if not re.fullmatch(r'[A-Z]{2}_[A-Z0-9]+\.[1-9][0-9]*', accession):
                        raise ValueError(f'{name}: expected an exact RefSeq accession.version header')
                    # The prefix filter only skips catalogue records which cannot be in
                    # any inventory. Exact accession.version is always the join key.
                    prefixes.add(accession.partition('_')[0])
                    batch.append((name, col, seq, accession))
                    count += 1
                    if len(batch) == 50000:
                        _insert_records(db, batch)
                        batch.clear()
                _insert_records(db, batch)
            report['chunks'][name] = {'records': count, 'inventory_sha256': sha256(inventory)}
            print(json.dumps({'event': 'inventory', 'chunk': name, 'records': count}), flush=True)
        opener = gzip.open if catalogue.suffix == '.gz' else open
        with opener(catalogue, 'rt', encoding='utf-8') as stream:
            batch = []
            for line_number, line in enumerate(stream, 1):
                if line_number % 1000000 == 0:
                    print(json.dumps({'event': 'catalogue_scan', 'rows_scanned': line_number}), flush=True)
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 6:
                    raise ValueError(f'Catalogue row {line_number}: expected six columns')
                accession = fields[2]
                if accession.partition('_')[0] not in prefixes:
                    continue
                batch.append((number(fields[0], positive=True), accession))
                if len(batch) == 50000:
                    db.executemany('UPDATE accession SET taxid=? WHERE name=?', batch)
                    db.commit()
                    batch.clear()
            db.executemany('UPDATE accession SET taxid=? WHERE name=?', batch)
            db.commit()
        missing = db.execute('SELECT COUNT(*) FROM accession WHERE taxid IS NULL').fetchone()[0]
        report['unique_accessions'] = db.execute('SELECT COUNT(*) FROM accession').fetchone()[0]
        report['missing_accessions'] = missing
        if missing:
            report['missing_examples'] = [row[0] for row in db.execute(
                'SELECT name FROM accession WHERE taxid IS NULL ORDER BY name LIMIT 20')]
            (output / 'coverage.json').write_text(json.dumps(report, indent=2) + '\n')
            raise ValueError(f'{missing} indexed accessions have no release-catalogue taxid; no sidecars published')
        # Publish only after every inventory has complete coverage.
        for name, info in report['chunks'].items():
            folder = output / name
            folder.mkdir()
            destination = folder / 'annotation.taxids.tsv'
            with destination.open('x', encoding='utf-8', newline='\n') as stream:
                stream.write('METAGRAPH_RECORD_TAXIDS_V1\n')
                count = 0
                for col, seq, accession, taxid in db.execute('''
                        SELECT r.col,r.seq,r.accession,a.taxid FROM record r
                        JOIN accession a ON a.name=r.accession
                        WHERE r.chunk=? ORDER BY r.col,r.seq''', (name,)):
                    stream.write(f'{col}\t{seq}\t{accession}\t{taxid}\n')
                    count += 1
            if count != info['records']:
                raise RuntimeError('Join lost indexed records')
            info['mapped_records'] = count
            info['sidecar_sha256'] = sha256(destination)
        report['complete'] = True
        report['chunks_checked'] = len(chunks)
        (output / 'coverage.json').write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps({'event': 'complete', 'chunks_checked': len(chunks),
                          'chunks_total': len(chunks), 'missing_accessions': 0}), flush=True)
    return report


def _insert_records(db, batch):
    db.executemany('INSERT INTO record VALUES (?,?,?,?)', batch)
    db.executemany('INSERT OR IGNORE INTO accession(name) VALUES (?)', ((r[3],) for r in batch))
    db.commit()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path, required=True)
    parser.add_argument('--chunk', action='append', required=True, metavar='NAME=HEADERS_TSV')
    parser.add_argument('--output', type=Path, required=True, help='New output directory; parent must exist')
    args = parser.parse_args()
    chunks = {}
    for item in args.chunk:
        name, separator, path = item.partition('=')
        if not separator or name in chunks:
            parser.error('Each chunk needs a unique NAME=HEADERS_TSV')
        chunks[name] = path
    prepare(args.catalogue, chunks, args.output)


if __name__ == '__main__':
    main()
