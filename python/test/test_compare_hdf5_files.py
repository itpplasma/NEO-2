import h5py
import pytest

from neo2_util import compare_hdf5_files


@pytest.mark.parametrize("group", ["fulltransp", "neo2_multispecies_out"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("value, expected", [(2.0, True), (2.001, False)])
def test_fulltransp_renamed_dataset(tmp_path, group, reverse, value, expected):
    reference = tmp_path / "reference.h5"
    other = tmp_path / "other.h5"
    names = ["avnabpsi", "av_nabla_stor"]
    if reverse:
        names.reverse()
    for filename, name, data in zip([reference, other], names, [2.0, value]):
        with h5py.File(filename, "w") as file:
            file.create_dataset(group + "/" + name, data=data)
    assert compare_hdf5_files(reference, other, 1e-8, [], [], True) == [expected, True]


@pytest.mark.parametrize("file_with_both_names", ["reference", "other"])
def test_file_with_old_and_new_name_is_rejected(tmp_path, file_with_both_names):
    # If one file holds both avnabpsi and av_nabla_stor, it is unclear which
    # of the two to compare with the single dataset of the other file, so the
    # comparison must refuse instead of picking one.
    paths = {"reference": tmp_path / "reference.h5", "other": tmp_path / "other.h5"}
    for role, path in paths.items():
        with h5py.File(path, "w") as file:
            file.create_dataset("fulltransp/av_nabla_stor", data=2.0)
            if role == file_with_both_names:
                file.create_dataset("fulltransp/avnabpsi", data=2.0)
    with pytest.raises(ValueError, match="Ambiguous dataset names"):
        compare_hdf5_files(paths["reference"], paths["other"], 1e-8, [], [], True)


def test_renamed_dataset_reference_accuracy_and_filter(tmp_path):
    reference = tmp_path / "reference.h5"
    other = tmp_path / "other.h5"
    with h5py.File(reference, "w") as file:
        dataset = file.create_dataset("fulltransp/avnabpsi", data=2.0)
        dataset.attrs["accuracy"] = 0.01
    with h5py.File(other, "w") as file:
        file.create_dataset("fulltransp/av_nabla_stor", data=2.001)
    assert compare_hdf5_files(reference, other, 1e-8, [], [], True) == [True, True]
    with h5py.File(reference, "r+") as file:
        file["fulltransp/avnabpsi"].attrs["accuracy"] = 1e-8
    assert compare_hdf5_files(reference, other, 1e-8, ["fulltransp", "avnabpsi"], [], True) == [False, True]
    assert compare_hdf5_files(reference, other, 1e-8, [], ["avnabpsi"], True) == [True, True]


def test_rename_is_limited_to_fulltransp(tmp_path):
    reference = tmp_path / "reference.h5"
    other = tmp_path / "other.h5"
    for filename, name in [(reference, "avnabpsi"), (other, "av_nabla_stor")]:
        with h5py.File(filename, "w") as file:
            file.create_dataset("unrelated/" + name, data=2.0)
    assert compare_hdf5_files(reference, other, 1e-8, [], [], True)[1] is False
