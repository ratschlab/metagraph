"""Case must remain meaningful through construction, annotation and lookup."""
import json
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import unittest

from base import METAGRAPH


@unittest.skipUnless('DNA_CASE_SENSITIVE' in os.path.realpath(METAGRAPH),
                     'requires the case-sensitive executable')
class TestCaseSensitive(unittest.TestCase):
    def test_complements_preserve_case(self):
        complement = str.maketrans('ACGTNacgt', 'TGCANtgca')
        sequences = ['ACG', 'acg', 'AcG', 'ANT', 'aNt']
        reverse = lambda sequence: sequence.translate(complement)[::-1]
        queries = sorted(set(sequences + [reverse(s) for s in sequences]))
        for backend, extension in [('hashfast', '.hashfastdbg'), ('succinct', '.dbg'),
                                   ('sshash', '.sshashdbg')]:
            for mode in ['basic', 'canonical', 'primary']:
                with self.subTest(backend=backend, mode=mode), TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    fasta = root / 'input.fa'
                    fasta.write_text(''.join(f'>{s}\n{s}\n' for s in sequences))
                    query = root / 'query.fa'
                    query.write_text(''.join(f'>{s}\n{s}\n' for s in queries))
                    graph = root / 'graph'
                    annotation = root / 'annotation'
                    def run(*args):
                        result = subprocess.run([METAGRAPH, *map(str, args)],
                                                capture_output=True, text=True)
                        self.assertEqual(result.returncode, 0, result.stderr)
                        return result.stdout
                    run('build', '--graph', backend, '--mode', mode, '-k', '3', '-p', '1',
                        '-o', graph, fasta)
                    graph = str(graph) + extension
                    run('annotate', '-i', graph, '--anno-header', '-o', annotation, fasta)
                    output = run('query', '-i', graph, '-a', str(annotation) + '.column.annodbg',
                                 '--query-mode', 'labels', '--min-kmers-fraction-label', '1',
                                 '--json', query)
                    results = {}
                    decoder = json.JSONDecoder()
                    while output.strip():
                        output = output.lstrip()
                        item, end = decoder.raw_decode(output)
                        output = output[end:]
                        results[item['seq_description']] = {r['sample'] for r in item['results']}
                    self.assertEqual(set(results), set(queries))
                    for q in queries:
                        expected = {s for s in sequences
                                    if q == s or (mode != 'basic' and q == reverse(s))}
                        self.assertEqual(results[q], expected, q)
