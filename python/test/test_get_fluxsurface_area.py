import pytest
import h5py
import numpy as np

@pytest.fixture
def neo2_output():
    filename = "neo2_ouput.h5"
    clean_up(filename)
    n_surf = 10
    content = {}
    content["boozer_s"] = np.linspace(0, 1, n_surf)
    content["avnabstor"] = np.linspace(0.1, 1, n_surf)
    content["avbhat2"] = np.linspace(1, 0.91, n_surf)
    content["Bref"] = np.linspace(0.1, 1, n_surf)
    content["psi_pr_hat"] = -np.linspace(0.1, 1, n_surf)
    content["bcovar_tht"] = -np.linspace(1, 0.1, n_surf)
    content["bcovar_phi"] = -np.linspace(0.1, 1, n_surf)
    content["aiota"] = np.linspace(1, 0.1, n_surf)
    with h5py.File(filename, "w") as file:
        for name in content.keys():
            dataset = file.create_dataset(name,(n_surf,), dtype="f")
            dataset[:] = content[name]
    yield filename, content
    clean_up(filename)

@pytest.fixture
def neo2_output_single_surface():
    filename = "neo2_ouput.h5"
    clean_up(filename)
    content = {}
    content["boozer_s"] = 0.5
    content["avnabstor"] = 0.5
    content["avbhat2"] = 0.5 
    content["Bref"] = 0.5
    content["psi_pr_hat"] = -0.5
    content["bcovar_tht"] = -0.5
    content["bcovar_phi"] = -0.5
    content["aiota"] = 0.5
    with h5py.File(filename, "w") as file:
        for name in content.keys():
            dataset = file.create_dataset(name, data=content[name], dtype="f")
    yield filename, content
    clean_up(filename)


def clean_up(filename):
    import os
    if os.path.exists(filename):
        os.remove(filename)


def test_read_neo2_output(neo2_output):
    filename, content = neo2_output
    with h5py.File(filename, "r") as file:
        for name in content.keys():
            assert np.allclose(content[name], file[name][:])


def test_call_get_fluxsurface_area(neo2_output):
    from neo2_ql import get_fluxsurface_area
    filename, _ = neo2_output
    _, _ = get_fluxsurface_area(filename)


def test_call_get_fluxsurface_area_single_surface(neo2_output_single_surface):
    from neo2_ql import get_fluxsurface_area
    filename, _ = neo2_output_single_surface
    _, _ = get_fluxsurface_area(filename)


def test_get_average_nable_stor(neo2_output):
    from neo2_ql import get_average_nabla_stor
    filename, content = neo2_output
    average_nabla_stor, stor = get_average_nabla_stor(filename)
    assert np.allclose(average_nabla_stor, content["avnabstor"])
    assert np.allclose(stor, content["boozer_s"])

    
