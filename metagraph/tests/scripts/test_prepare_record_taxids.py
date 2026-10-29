import gzip
import importlib.util
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / 'scripts/prepare_record_taxids.py'
spec = importlib.util.spec_from_file_location('prepare_record_taxids', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class RecordTaxidPreparation(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def inventory(self, name, rows):
        path = self.root / (name + '.tsv')
        path.write_text('METAGRAPH_SEQUENCE_HEADERS_V1\n' + ''.join(rows))
        return path

    def catalogue(self, rows):
        path = self.root / 'RefSeq-release237.catalog.gz'
        with gzip.open(path, 'wt') as stream:
            for accession, taxid in rows:
                stream.write(f'{taxid}\torganism\t{accession}\tbacteria|complete\tna\t1000\n')
        return path

    def test_complete_join_preserves_every_record_and_duplicate_accession(self):
        chunks = {
            'first': self.inventory('first', ['0\t0\tNC_000001.1\n', '0\t1\tNZ_ABC001.2\n']),
            'second': self.inventory('second', ['0\t0\tNC_000001.1\n']),
        }
        catalog = self.catalogue([('WP_999999.1', 2), ('NC_000001.1', 100), ('NZ_ABC001.2', 200)])
        output = self.root / 'output'
        report = module.prepare(catalog, chunks, output)
        self.assertTrue(report['complete'])
        self.assertEqual((2, 2, 2, 0), (report['chunks_total'], report['chunks_checked'],
                                        report['unique_accessions'], report['missing_accessions']))
        self.assertEqual('METAGRAPH_RECORD_TAXIDS_V1\n0\t0\tNC_000001.1\t100\n0\t1\tNZ_ABC001.2\t200\n',
                         (output / 'first/annotation.taxids.tsv').read_text())
        self.assertEqual(1, report['chunks']['second']['mapped_records'])
        self.assertEqual(module.sha256(output / 'second/annotation.taxids.tsv'),
                         report['chunks']['second']['sidecar_sha256'])
        with self.assertRaises(FileExistsError):
            module.prepare(catalog, chunks, output)

    def test_missing_version_does_not_publish_partial_mapping(self):
        headers = self.inventory('first', ['0\t0\tNC_000001.2\n'])
        catalog = self.catalogue([('NC_000001.1', 100)])
        output = self.root / 'missing'
        with self.assertRaisesRegex(ValueError, 'no release-catalogue taxid'):
            module.prepare(catalog, {'first': headers}, output)
        self.assertEqual([], list(output.rglob('annotation.taxids.tsv')))
        report = json.loads((output / 'coverage.json').read_text())
        self.assertFalse(report['complete'])
        self.assertEqual(['NC_000001.2'], report['missing_examples'])

    def test_conflicting_taxids_and_invalid_taxid_are_refused(self):
        headers = self.inventory('first', ['0\t0\tNC_000001.1\n'])
        for i, rows in enumerate([
                [('NC_000001.1', 100), ('NC_000001.1', 200)], [('NC_000001.1', 0)]]):
            with self.assertRaises((sqlite3.IntegrityError, ValueError)):
                module.prepare(self.catalogue(rows), {'first': headers}, self.root / f'bad{i}')
            self.assertEqual([], list((self.root / f'bad{i}').rglob('annotation.taxids.tsv')))

    def test_missing_record_index_is_refused(self):
        headers = self.inventory('first', ['0\t0\tNC_000001.1\n', '0\t2\tNC_000002.1\n'])
        catalog = self.catalogue([('NC_000001.1', 100), ('NC_000002.1', 200)])
        with self.assertRaisesRegex(ValueError, 'missing record ids'):
            module.prepare(catalog, {'first': headers}, self.root / 'gap')


if __name__ == '__main__':
    unittest.main()
