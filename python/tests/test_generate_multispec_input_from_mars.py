#%% Standard modules
import numpy as np
import os
import h5py

# Modules to test
from neo2_mars import get_species_cgs_from_mars
from neo2_mars import generate_multispec_input_from_mars

test_mars_dir = '/proj/plasma/DATA/DEMO/MARS/MARSQ_INPUTS_KNTV21_NEO2profs_RUN/'
test_outpu_dir = '/tmp/'

def test_get_species_cgs_from_mars():
    Ze, Zi, me, mi = get_species_cgs_from_mars(test_mars_dir)
    assert Ze == -1
    assert Zi == 1
    assert np.isclose(me, 9.10938356e-28, rtol=1e-4, atol=0)
    assert np.isclose(mi, 3.343583719e-24, rtol=1e-3, atol=0) # off due MARS using not atomic unit, but proton mass

def test_generate_multispec_input_from_mars():
    output_file = os.path.join(test_outpu_dir, 'multi_spec.in')
    generate_multispec_input_from_mars(test_mars_dir, output_file, 
                                       number_of_surfaces=10, bounds=[0.0, 1.0])
    print_hdf5_structure(output_file)

def print_hdf5_structure(file_name):
    def print_attrs(name, obj):
        print(f"\n{name}")
        for key, val in obj.attrs.items():
            print(f"    Attribute: {key} => {val}")

    def print_dataset(name, dataset):
        print(f"\nDataset: {name}")
        print(f"    Shape: {dataset.shape}")
        print(f"    Data: {dataset[...]}")
    
    with h5py.File(file_name, 'r') as file:
        print(f"HDF5 file: {file_name}")
        file.visititems(lambda name, obj: print_attrs(name, obj) if isinstance(obj, h5py.Group) else print_dataset(name, obj))


if __name__ == '__main__':
    test_get_species_cgs_from_mars()
    test_generate_multispec_input_from_mars()
    print('All tests passed.')