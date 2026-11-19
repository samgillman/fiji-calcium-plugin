import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
import build_release


class ReleaseTests(unittest.TestCase):
    def test_package_is_installable_and_reproducible(self):
        build_release.build()
        archive = build_release.ROOT / 'dist/Calcium_Imaging_Fiji_Plugin.zip'
        first = archive.read_bytes()
        build_release.build()
        self.assertEqual(first, archive.read_bytes())
        with zipfile.ZipFile(archive) as z:
            self.assertIsNone(z.testzip())
            self.assertIn(build_release.PLUGIN.as_posix(), z.namelist())
            self.assertEqual(z.read(build_release.PLUGIN.as_posix()), (build_release.ROOT / build_release.PLUGIN).read_bytes())
            self.assertFalse(any(name.startswith('/') or '..' in Path(name).parts for name in z.namelist()))
        for line in (archive.parent / 'SHA256SUMS.txt').read_text().splitlines():
            digest, name = line.split('  ')
            self.assertEqual(digest, hashlib.sha256((archive.parent / name).read_bytes()).hexdigest())

    def test_mismatched_version_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'VERSION').write_text('2026.11.20.1\n')
            script = root / build_release.PLUGIN
            script.parent.mkdir(parents=True)
            script.write_text('static final String VERSION = "2026.11.19.2"')
            with patch.object(build_release, 'ROOT', root):
                with self.assertRaises(ValueError):
                    build_release.validate()


if __name__ == '__main__':
    unittest.main()
