"""Independent malformed-output oracles for the dependency golden guard."""

import importlib.util
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

script = Path(__file__).resolve().parents[1] / "ShellScripts/check_golden_hdf5.py"
spec = importlib.util.spec_from_file_location("golden_guard", script)
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


class GoldenHdf5GuardTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        root = Path(directory.name)
        self.reference = root / "reference.h5"
        self.current = root / "current.h5"
        for path in (self.reference, self.current):
            with h5py.File(path, "w") as file:
                file["transport/coefficient"] = [1.0, 2.0]
                file["response"] = np.array([1.0 + 2.0j])

    def test_identical_finite_outputs_pass(self):
        guard.check(self.reference, self.current, 1e-14)

    def test_missing_or_extra_dataset_fails(self):
        with h5py.File(self.current, "a") as file:
            del file["transport/coefficient"]
        with self.assertRaises(ValueError):
            guard.check(self.reference, self.current, 1e-14)
        with self.assertRaises(ValueError):
            guard.check(self.current, self.reference, 1e-14)

    def test_shape_change_fails(self):
        with h5py.File(self.current, "a") as file:
            del file["transport/coefficient"]
            file["transport/coefficient"] = [[1.0, 2.0]]
        with self.assertRaises(ValueError):
            guard.check(self.reference, self.current, 1e-14)

    def test_nan_and_infinity_fail_in_either_output(self):
        for path in (self.reference, self.current):
            for value in (np.nan, np.inf, -np.inf):
                with h5py.File(path, "a") as file:
                    file["transport/coefficient"][0] = value
                with self.assertRaises(ValueError):
                    guard.check(self.reference, self.current, 1e-14)
            with h5py.File(path, "a") as file:
                file["transport/coefficient"][0] = 1.0

    def test_changed_complex_response_fails(self):
        with h5py.File(self.current, "a") as file:
            file["response"][0] = 1.0 + 3.0j
        with self.assertRaises(ValueError):
            guard.check(self.reference, self.current, 1e-14)


if __name__ == "__main__":
    unittest.main()
