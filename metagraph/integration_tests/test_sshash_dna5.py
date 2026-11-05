"""DNA5 self-complements retain their labels through SSHash indexing."""
import itertools
import json
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import unittest

from base import METAGRAPH


@unittest.skipUnless('DNA5' in os.path.realpath(METAGRAPH), 'requires the DNA5 executable')
class TestSSHashDNA5(unittest.TestCase):
    def test_all_three_mers(self):
        words = [''.join(w) for w in itertools.product('ACGTN', repeat=3)]
        reverse = lambda s: s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
        representatives = sorted({min(s, reverse(s)) for s in words})
        for mode in ['basic', 'canonical', 'primary']:
            with self.subTest(mode=mode), TemporaryDirectory() as tmp:
                root = Path(tmp)
                source, queries = root / 'input.fa', root / 'queries.fa'
                source.write_text(''.join(f'>{s}\n{s}\n' for s in representatives))
                queries.write_text(''.join(f'>{s}\n{s}\n' for s in words))
                graph, annotation = root / 'graph', root / 'annotation'
                def run(*args):
                    result = subprocess.run([METAGRAPH, *map(str, args)], capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    return result.stdout
                run('build', '--graph', 'sshash', '--mode', mode, '-k', '3', '-p', '4', '-o', graph, source)
                graph = str(graph) + '.sshashdbg'
                run('annotate', '-i', graph, '--anno-header', '-o', annotation, source)
                output = run('query', '-i', graph, '-a', str(annotation) + '.column.annodbg',
                             '--query-mode', 'labels', '--min-kmers-fraction-label', '1', '--json', queries)
                results, decoder = {}, json.JSONDecoder()
                while output.strip():
                    output = output.lstrip()
                    item, end = decoder.raw_decode(output)
                    output = output[end:]
                    results[item['seq_description']] = {r['sample'] for r in item['results']}
                self.assertEqual(set(results), set(words))
                for word in words:
                    expected = {word} if word in representatives else set()
                    if mode != 'basic':
                        expected = {min(word, reverse(word))}
                    self.assertEqual(results[word], expected, word)
