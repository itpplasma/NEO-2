"""Behavioral checks for cached package bytes and fresh APT metadata."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


HELPER = os.environ.get(
    "APT_CACHE_HELPER", str(Path(__file__).with_name("verify_apt_cache.py"))
)
PACKAGE = "libmpfr6_4.2.1-1build1_amd64.deb"
PACKAGE_2 = "libsuitesparseconfig7_1%3a7.6.1+dfsg-1build1_amd64.deb"
BYTES_A = b"!<arch>\nindependent known package fixture A\n"
BYTES_B = b"!<arch>\nindependent known package fixture B\n"
SHA_A = "52d7dc2c7e4ce072f8f627232e90f72b01c221f4820cf07ff2675e905ad8f4c7"
SHA_B = "5f60f4b83960d56041ca31b0317467358534be8b6b71b2d41442b80d336aa70c"


def package_line(name=PACKAGE, size=44, digest=SHA_A):
    return (
        "'https://archive.ubuntu.com/ubuntu/pool/main/p/package.deb' "
        f"{name} {size} SHA256:{digest}\n"
    )


class VerifyAptCacheTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.manifest = self.root / "manifest.txt"
        self.archives = self.root / "archives"
        self.archives.mkdir()

    def run_helper(self, text):
        self.manifest.write_text(text, encoding="utf-8")
        return subprocess.run(
            [sys.executable, HELPER, str(self.manifest), str(self.archives)],
            capture_output=True, text=True, timeout=10, check=False,
        )

    def assert_success(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertRegex(result.stdout.strip(), r"^cache-key=[0-9a-f]{64}$")

    def cached(self, name=PACKAGE, data=BYTES_A):
        path = self.archives / name
        path.write_bytes(data)
        return path

    def test_cold_cache_is_successful(self):
        self.assert_success(self.run_helper(package_line()))
        self.assertEqual(list(self.archives.glob("*.deb")), [])

    def test_matching_warm_cache_preserves_exact_bytes(self):
        cached = self.cached()
        self.assert_success(self.run_helper(package_line()))
        self.assertEqual(cached.read_bytes(), BYTES_A)

    def test_same_size_corruption_is_removed(self):
        corrupted = self.cached(data=BYTES_B)
        self.assertEqual(len(BYTES_A), len(BYTES_B))
        self.assert_success(self.run_helper(package_line()))
        self.assertFalse(corrupted.exists())

    def test_truncated_package_is_removed(self):
        corrupted = self.cached(data=BYTES_A[:-1])
        self.assert_success(self.run_helper(package_line()))
        self.assertFalse(corrupted.exists())

    def test_directory_named_like_package_is_discarded(self):
        invalid = self.archives / PACKAGE
        invalid.mkdir()
        (invalid / "unexpected-child").write_bytes(BYTES_B)
        partial = self.archives / "partial"
        partial.mkdir()
        (partial / "keep").write_bytes(BYTES_A)
        self.assert_success(self.run_helper(package_line()))
        self.assertFalse(invalid.exists())
        self.assertEqual((partial / "keep").read_bytes(), BYTES_A)

    @unittest.skipIf(hasattr(os, "geteuid") and os.geteuid() == 0,
                     "root can read mode-000 files")
    def test_unreadable_cached_package_is_discarded(self):
        invalid = self.cached()
        invalid.chmod(0)
        self.assert_success(self.run_helper(package_line()))
        self.assertFalse(invalid.exists())

    def test_updated_manifest_rejects_previously_valid_package(self):
        outdated = self.cached()
        self.assert_success(self.run_helper(package_line(digest=SHA_B)))
        self.assertFalse(outdated.exists())

    def test_partial_hit_keeps_valid_and_removes_unlisted_archives(self):
        valid = self.cached()
        outdated = self.cached("libmpfr6_old_amd64.deb", BYTES_B)
        malformed = self.cached("bad name.deb", BYTES_B)
        manifest = package_line() + package_line(PACKAGE_2, digest=SHA_B)
        self.assert_success(self.run_helper(manifest))
        self.assertEqual(valid.read_bytes(), BYTES_A)
        self.assertFalse(outdated.exists())
        self.assertFalse(malformed.exists())
        self.assertFalse((self.archives / PACKAGE_2).exists())

    def test_symlink_is_removed_without_modifying_its_target(self):
        outside = self.root / "outside.deb"
        outside.write_bytes(BYTES_A)
        cached = self.archives / PACKAGE
        cached.symlink_to(outside)
        self.assert_success(self.run_helper(package_line()))
        self.assertFalse(cached.is_symlink())
        self.assertEqual(outside.read_bytes(), BYTES_A)
        other_link = self.archives / "unexpected-link"
        other_link.symlink_to(outside)
        self.assert_success(self.run_helper(package_line()))
        self.assertFalse(other_link.is_symlink())
        self.assertEqual(outside.read_bytes(), BYTES_A)

    def test_apt_headers_and_epoch_encoded_filename_are_allowed(self):
        cached = self.cached(PACKAGE_2, BYTES_B)
        manifest = (
            "Reading package lists...\nBuilding dependency tree...\n"
            "The following NEW packages will be installed:\n  libsuitesparseconfig7\n"
            "0 upgraded, 1 newly installed, 0 to remove and 0 not upgraded.\n"
            + package_line(PACKAGE_2, digest=SHA_B)
        )
        self.assert_success(self.run_helper(manifest))
        self.assertEqual(cached.read_bytes(), BYTES_B)

    def test_already_installed_empty_manifest_removes_unneeded_cache(self):
        cached = self.cached()
        manifest = (
            "Reading package lists...\nlibmpfr6 is already the newest version.\n"
            "0 upgraded, 0 newly installed, 0 to remove and 0 not upgraded.\n"
        )
        self.assert_success(self.run_helper(manifest))
        self.assertFalse(cached.exists())

    def test_malformed_trusted_package_lines_fail_closed(self):
        invalid_lines = [
            package_line(size=-1), package_line(size="unknown"),
            package_line(digest="0" * 63), package_line(digest="g" * 64),
            package_line(digest=SHA_A.upper()),
            package_line().replace("SHA256:", "SHA1:"),
            package_line().replace("SHA256:", ""),
            package_line().rstrip() + " unexpected\n",
            "'https://archive.ubuntu.com/truncated.deb\n",
        ]
        for line in invalid_lines:
            with self.subTest(line=line):
                unlisted = self.cached("unlisted_1_amd64.deb", BYTES_B)
                result = self.run_helper(package_line() + line)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(unlisted.read_bytes(), BYTES_B)

    def test_manifest_filenames_cannot_escape_archives(self):
        for name in ("../outside.deb", "/absolute.deb", "nested/package.deb"):
            with self.subTest(name=name):
                result = self.run_helper(package_line(name=name))
                self.assertNotEqual(result.returncode, 0)

    def test_conflicting_manifest_entries_fail_closed(self):
        result = self.run_helper(package_line() + package_line(digest=SHA_B))
        self.assertNotEqual(result.returncode, 0)

    def test_cache_key_follows_plan_independent_of_order_and_cached_bytes(self):
        first_line = package_line()
        second_line = package_line(PACKAGE_2, digest=SHA_B)
        first = self.run_helper(first_line + second_line)
        self.assert_success(first)
        self.cached(data=BYTES_B)
        reordered = self.run_helper("Reading package lists...\n" + second_line + first_line)
        self.assert_success(reordered)
        self.assertEqual(first.stdout, reordered.stdout)
        changed = self.run_helper(package_line(digest=SHA_B) + second_line)
        self.assert_success(changed)
        self.assertNotEqual(first.stdout, changed.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
