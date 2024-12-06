import pytest
import h5py
import numpy as np

@pytest.fixture
def neo2_output():
    filename = "neo2_ouput.h5"
    clean_up(filename)
    n_surf = 10
    content = {}
    content["avnabpsi"] = np.linspace(0.1,1,n_surf)
    with h5py.File(filename, "w") as file:
        average_nabla_stor = file.create_dataset("avnabpsi",(n_surf,), dtype="f")
        average_nabla_stor[:] = content["avnabpsi"]
    yield filename, content
    clean_up(filename)

def clean_up(filename):
    import os
    if os.path.exists(filename):
        os.remove(filename)

def test_read_neo2_output(neo2_output):
    filename , content = neo2_output
    with h5py.File(filename, "r") as file:
        assert np.allclose(content["avnabpsi"], file["avnabpsi"][:])
