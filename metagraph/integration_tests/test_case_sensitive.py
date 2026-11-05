"""Case encodes strand throughout graph construction, lookup and extraction."""
import gzip
import itertools
import json
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import unittest

from base import METAGRAPH


COMPLEMENT = str.maketrans("ACGTNacgt", "tgcaNTGCA")


def reverse(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def kmers(sequences):
    return {s[i:i + 3] for s in sequences for i in range(len(s) - 2)}


@unittest.skipUnless('DNA_CASE_SENSITIVE' in os.path.realpath(METAGRAPH),
                     'requires the case-sensitive executable')
class TestCaseSensitive(unittest.TestCase):
    def run_cli(self, *args):
        result = subprocess.run([METAGRAPH, *map(str, args)], capture_output=True,
                                text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout

    def check_queries(self, graph, root, fasta, query, sequences, both):
        annotation = root / 'annotation'
        self.run_cli('annotate', '-i', graph, '--anno-header', '-o', annotation,
                     *(['--fwd-and-reverse'] if both == 'input' else []), fasta)
        output = self.run_cli('query', '-i', graph, '-a', str(annotation) + '.column.annodbg',
                              '--query-mode', 'labels', '--min-kmers-fraction-label', '1',
                              '--json', query)
        results = {}
        decoder = json.JSONDecoder()
        while output.strip():
            item, end = decoder.raw_decode(output.lstrip())
            output = output.lstrip()[end:]
            results[item['seq_description']] = {r['sample'] for r in item['results']}
        self.assertEqual(len(results), 9 ** 3)
        for q, labels in results.items():
            expected = {s for s in sequences
                        if q in s or (both and reverse(q) in s)}
            self.assertEqual(labels, expected, q)

    def test_strand_complements_and_extraction(self):
        sequences = ['ACG', 'acg', 'AcG', 'ANT', 'aNt', 'NNN',
                     'ACGTAC', 'acgtac', 'AcGNta']
        queries = [''.join(s) for s in itertools.product('ACGTNacgt', repeat=3)]
        for backend, extension in [('hashfast', '.hashfastdbg'), ('succinct', '.dbg')]:
            for mode in ['basic', 'canonical', 'primary']:
                with self.subTest(backend=backend, mode=mode), TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    fasta, query = root / 'input.fa', root / 'query.fa'
                    fasta.write_text(''.join(f'>{s}\n{s}\n' for s in sequences))
                    query.write_text(''.join(f'>{s}\n{s}\n' for s in queries))
                    if mode == 'primary':
                        # Primary mode expects one representative per strand pair.
                        self.run_cli('build', '--graph', backend, '--mode', 'canonical',
                                     '-k', '3', '-p', '1', '-o', root / 'canonical', fasta)
                        self.run_cli('transform', '--to-fasta', '--primary-kmers',
                                     '-o', root / 'primary', str(root / 'canonical') + extension)
                        build_input = root / 'primary.fasta.gz'
                    else:
                        build_input = fasta
                    self.run_cli('build', '--graph', backend, '--mode', mode, '-k', '3',
                                 '-p', '1', '-o', root / 'graph', build_input)
                    graph = str(root / 'graph') + extension
                    self.check_queries(graph, root, fasta, query, sequences, mode != 'basic')
                    normalize = (lambda s: min(s, reverse(s))) if mode != 'basic' else (lambda s: s)
                    expected = {normalize(s) for s in kmers(sequences)}
                    # Check both contig and unitig extraction, including strand deduplication.
                    for unitigs in [False, True]:
                        prefix = root / ('unitigs' if unitigs else 'contigs')
                        self.run_cli('transform', '--to-fasta',
                                     *(['--unitigs'] if unitigs else []),
                                     *(['--primary-kmers'] if mode != 'basic' else []),
                                     '-o', prefix, graph)
                        exported = Path(str(prefix) + '.fasta.gz')
                        with gzip.open(exported, 'rt') as stream:
                            records = stream.read().split('>')[1:]
                        extracted = [''.join(record.splitlines()[1:]) for record in records]
                        self.assertEqual({normalize(s) for s in kmers(extracted)}, expected)
                        # Rebuild exported sequences and annotate using the original input.
                        rebuilt = str(prefix) + '-rebuilt'
                        self.run_cli('build', '--graph', backend, '--mode', mode, '-k', '3',
                                     '-p', '1', '-o', rebuilt, exported)
                        self.check_queries(rebuilt + extension, root, fasta, query,
                                           sequences, mode != 'basic')

    def test_reverse_input_in_basic_graph(self):
        sequences = ['ACG', 'AcGNta']
        for backend, extension in [('hashfast', '.hashfastdbg'), ('succinct', '.dbg')]:
            with self.subTest(backend=backend), TemporaryDirectory() as tmp:
                root = Path(tmp)
                fasta, query = root / 'input.fa', root / 'query.fa'
                fasta.write_text(''.join(f'>{s}\n{s}\n' for s in sequences))
                query.write_text(''.join(f'>{s}\n{s}\n' for s in
                                        map(''.join, itertools.product('ACGTNacgt', repeat=3))))
                self.run_cli('build', '--graph', backend, '--mode', 'basic',
                             '--fwd-and-reverse', '-k', '3', '-p', '1',
                             '-o', root / 'graph', fasta)
                self.check_queries(str(root / 'graph') + extension, root, fasta, query,
                                   sequences, 'input')
