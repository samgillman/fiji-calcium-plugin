import base64
import hashlib
import io
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
import embed_analysis


class AnalysisBundleTests(unittest.TestCase):
    def test_upstream_matches_recorded_provenance(self):
        root = embed_analysis.ROOT / 'analysis'
        provenance = dict(line.split('\t') for line in (root / 'provenance.tsv').read_text().splitlines())
        for name in ['utils.R', 'responders.R', 'response_protocols.R', 'LICENSE']:
            self.assertEqual(hashlib.sha256((root / 'upstream' / name).read_bytes()).hexdigest(), provenance[name + '_sha256'])

    def test_host_metadata_differences_pass_but_changed_sources_fail(self):
        source = embed_analysis.PLUGIN.read_text(encoding='utf-8')
        def payload(files):
            output = io.BytesIO()
            with zipfile.ZipFile(output, 'w') as archive:
                for name, data in files.items():
                    entry = zipfile.ZipInfo(name)
                    entry.create_system = 0  # Windows host, deliberately different metadata/compression.
                    archive.writestr(entry, data)
            data = output.getvalue()
            changed = re.sub(r'HASH = "[a-f0-9]{64}"', 'HASH = "' + hashlib.sha256(data).hexdigest() + '"', source)
            return re.sub(r'DATA = "[A-Za-z0-9+/=]+"', 'DATA = "' + base64.b64encode(data).decode() + '"', changed)
        with tempfile.TemporaryDirectory() as directory:
            plugin = Path(directory) / 'plugin.groovy'
            files = embed_analysis.sources()
            plugin.write_text(payload(files), encoding='utf-8')
            with patch.object(embed_analysis, 'PLUGIN', plugin):
                embed_analysis.sync(check=True)
                files['runner.R'] += b'\n# Unexpected change\n'
                plugin.write_text(payload(files), encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'embed_analysis.py'):
                    embed_analysis.sync(check=True)
