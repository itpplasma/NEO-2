# %%
from .get_integral_ntv import get_volume_element

def main():
    return

def get_average_nabla_stor(output_filename):
    import h5py
    with h5py.File(output_filename) as file:
        stor = file["boozer_s"][:]
        average_nabla_stor = file["avnabpsi"][:]
    return stor, average_nabla_stor

if __name__=="__main__":
    main()