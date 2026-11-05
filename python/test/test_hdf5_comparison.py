import h5py
import numpy as np
import pytest

from neo2_util import compare_hdf5_files


@pytest.fixture
def files(tmp_path):
    with h5py.File(tmp_path / "reference.h5", "w") as reference, \
            h5py.File(tmp_path / "other.h5", "w") as other:
        yield reference, other


def compare(files, delta=1e-3, whitelist=None, blacklist=None, verbose=False):
    for file in files:
        file.flush()
    return compare_hdf5_files(
        *(file.filename for file in files), delta,
        whitelist or [], blacklist or [], verbose
    )


@pytest.mark.parametrize("reference,other,equal", [
    ([1000., 0.001], [1000., 0.5], True),
    ([1000., 0.001], [1000., 2.], False),
    ([3. + 4.j, 0.j], [3. + 4.j, 0.004j], True),
    ([3. + 4.j, 0.j], [3. + 4.j, 0.006j], False),
    ([0., 0.], [0., 0.0005], True),
    ([0., 0.], [0., 0.002], False),
    ([], [], True),
    (1., 1.0005, True),
])
def test_dataset_max_tolerance(files, reference, other, equal):
    files[0]["nested/data"] = reference
    files[1]["nested/data"] = other
    assert compare(files) == [equal, True]


@pytest.mark.parametrize("accuracy,equal", [
    (0.1, False), (0.5, True), (0., False),
    (-1., True), (np.nan, True), (np.inf, True),
])
def test_absolute_accuracy_overrides_relative_tolerance(files, accuracy, equal):
    files[0]["data"] = [1000., 1.]
    files[0]["data"].attrs["accuracy"] = accuracy
    files[1]["data"] = [1000., 1.25]
    assert compare(files) == [equal, True]


@pytest.mark.parametrize("reference,other", [
    ([1., 2.], [[1., 2.]]),
    (1., [1.]),
    (np.empty((0, 2)), np.empty((0, 3))),
])
def test_shapes_must_match_without_broadcasting(files, reference, other):
    files[0]["data"] = reference
    files[1]["data"] = other
    assert compare(files) == [False, True]


@pytest.mark.parametrize("dataset_side", [0, 1])
def test_dataset_and_group_with_same_name_differ(files, dataset_side):
    files[dataset_side]["data"] = [1.]
    files[1 - dataset_side].create_group("data")
    assert compare(files) == [False, False]


@pytest.mark.parametrize("bad_value", [np.nan, np.inf, -np.inf, complex(1., np.nan)])
@pytest.mark.parametrize("bad_side", [0, 1])
def test_nonfinite_values_fail_in_either_file(files, bad_value, bad_side):
    files[bad_side]["nested/data"] = [bad_value]
    files[1 - bad_side]["nested/data"] = [1.]
    assert compare(files) == [False, True]


@pytest.mark.parametrize("bad_value", [np.nan, np.inf])
def test_matching_nonfinite_values_still_fail(files, bad_value):
    for file in files:
        file["data"] = [bad_value]
    assert compare(files) == [False, True]


@pytest.mark.parametrize("side", [0, 1])
@pytest.mark.parametrize("extra_value,values_equal", [(2., True), (np.nan, False)])
def test_non_overlapping_datasets_are_checked_for_nonfinite_values(files, side, extra_value, values_equal):
    for file in files:
        file["shared"] = [1.]
    files[side]["extra/nested/data"] = [extra_value]
    assert compare(files) == [values_equal, side == 1]


@pytest.mark.parametrize("reference,other", [
    (np.array([1], dtype=np.uint64), np.array([2], dtype=np.uint64)),
    (np.array([True]), np.array([False])),
    (np.array([1], dtype=np.int64), np.array([2], dtype=np.int64)),
    (np.array([b"first"]), np.array([b"second"])),
    (np.array(["first"], dtype=h5py.string_dtype()), np.array(["second"], dtype=h5py.string_dtype())),
    ([1.], [b"1"]),
])
def test_nonfloating_values_are_not_ignored(files, reference, other):
    files[0]["data"] = reference
    files[1]["data"] = other
    assert compare(files) == [False, True]


def test_numeric_storage_precision_need_not_match(files):
    files[0]["float"] = np.array([1., 2.], dtype=np.float64)
    files[1]["float"] = np.array([1., 2.], dtype=np.float32)
    files[0]["integer"] = np.array([1, 2], dtype=np.int64)
    files[1]["integer"] = np.array([1, 2], dtype=np.int32)
    assert compare(files) == [True, True]


def test_equal_nonfloating_values_pass(files):
    for file in files:
        file["unsigned"] = np.array([2**64 - 1], dtype=np.uint64)
        file["boolean"] = [True, False]
        file["fixed_string"] = [b"value"]
        file["variable_string"] = np.array(["value"], dtype=h5py.string_dtype())
    assert compare(files) == [True, True]


@pytest.mark.parametrize("verbose", [False, True])
def test_filters_apply_to_nested_and_non_overlapping_data(files, tmp_path, verbose):
    for file in files:
        file["keep/value"] = [1.]
        file["keep/ignore"] = [np.nan]
    files[1]["ignore/extra"] = [np.inf]
    blacklist = tmp_path / "blacklist.txt"
    blacklist.write_text("ignore\n")
    assert compare(files, blacklist=str(blacklist), verbose=verbose) == [True, True]
    assert compare(files, whitelist=["keep", "value"], verbose=verbose) == [True, True]


def test_filters_do_not_hide_missing_reference_keys(files):
    files[0]["ignored"] = [1.]
    assert compare(files, blacklist=["ignored"]) == [True, False]
