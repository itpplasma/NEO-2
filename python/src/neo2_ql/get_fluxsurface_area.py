# %%
from .get_integral_ntv import get_volume_element

def get_fluxsurface_area(neo2_ouput):
    dV_dstor, stor = get_volume_element(neo2_ouput)
    average_nabla_stor, _ = get_average_nabla_stor(neo2_ouput)
    area = dV_dstor * average_nabla_stor
    return area, stor

def get_average_nabla_stor(neo2_ouput):
    import h5py
    from neo2_util import get_hdf5dataset_value
    with h5py.File(neo2_ouput) as file:
        average_nabla_stor = get_hdf5dataset_value(file, "avnabpsi")
        stor = get_hdf5dataset_value(file, "boozer_s")
    return average_nabla_stor, stor

if __name__=="__main__":
    get_fluxsurface_area()