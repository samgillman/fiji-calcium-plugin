from pathlib import Path
import subprocess
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import release_status


class ReleaseStatusTests(unittest.TestCase):
    def test_unreleased_version_and_first_release(self):
        def lookup(repo, suffix):
            return None if suffix.startswith("tags/") else {"tag_name": "v2026.11.19.3"}
        self.assertTrue(release_status.needed("2026.11.19.4", "owner/repo", lookup))
        self.assertTrue(release_status.needed("2026.11.19.4", "owner/repo", lambda *_: None))
        with self.assertRaises(RuntimeError):
            release_status.needed("2026.11.19.2", "owner/repo", lookup)

    def test_existing_release_is_immutable_and_draft_blocks(self):
        self.assertFalse(release_status.needed("2026.11.19.4", "owner/repo", lambda *_: {"draft": False}))
        for existing in [{"draft": True}, {"prerelease": True}]:
            with self.assertRaises(RuntimeError):
                release_status.needed("2026.11.19.4", "owner/repo", lambda *_: existing)

    def test_only_real_404_means_absent(self):
        def response(code, stdout, stderr=""):
            return lambda *a, **k: subprocess.CompletedProcess(a, code, stdout, stderr)
        self.assertIsNone(release_status.release("owner/repo", "latest", response(1, '{"message":"Not Found","status":"404"}')))
        for body in ['{"message":"API rate limit exceeded","status":"403"}', '', '{"message":"Not Found","status":"500"}']:
            with self.assertRaises(RuntimeError):
                release_status.release("owner/repo", "latest", response(1, body, "request failed"))
        self.assertEqual(release_status.release("owner/repo", "latest", response(0, '{"tag_name":"v2026.11.19.4"}'))["tag_name"], "v2026.11.19.4")


if __name__ == "__main__":
    unittest.main()
