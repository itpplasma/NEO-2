"""
neo2_ql
-----------

This is the package containing tools for post-processing NEO-2-QL runs.
"""
from importlib import import_module

from .get_integral_ntv import get_integral_ntv_torque_neo2ql
from .load_profile_data import load_cgs_profiles_and_interp
from .load_profile_data import convert_units_from_norm_to_si, convert_units_from_si_to_cgs
from .load_profile_data import interp_grid, interp_profiles, interp_cubic
from .generate_multispec_input import generate_multispec_input
from .generate_multispec_input import write_multispec_to_hdf5
from .generate_multispec_input import get_coulomb_logarithm, get_kappa, derivative
from .generate_multispec_input import get_species_def_array
from .interpolate_profiles import interpolate_profiles_to_same_grid
from .neo2_output_omte import compute_neo2_er_from_transport_coefficients
from .neo2_output_omte import compute_neo2_omte_from_transport_coefficients
from .neo2_output_omte import compute_d31_reference_electron
from .neo2_output_omte import compute_omte_from_neo2_output
from .write_profiles import write_profiles_to_dat_files
from .get_fluxsurface_area import get_average_nabla_stor, get_fluxsurface_area
from .compute_omte import compute_omte_diamagnetic
from .compute_omte import compute_omte_force_balance, compute_omte_toroidal_rotation
from .compute_omte import compute_omte_toroidal_rotation_neo2_convention
from .compute_omte import compute_omte_neo2_single_ion_limit
from .compute_omte import compute_poloidal_rotation_neoclassical
from .compute_omte import compute_omte_neoclassical_poloidal
from .compute_omte import compute_omte_neoclassical_poloidal_auto_k
from .compute_omte import select_poloidal_rotation_coefficient
from .compute_omte import compute_boozer_metric
from .compute_omte import compute_boozer_metric_from_bc
from .compute_omte import compute_boozer_metric_from_rz_profile
from .compute_omte import compute_omte_toroidal_rotation_boozer
from .compute_omte import compute_omte_neoclassical_poloidal_boozer

_LAZY_EXPORTS = {
    'get_neo2_ql_input_profiles': ('.plot_neo2_ql_input_profiles', 'get_neo2_ql_input_profiles'),
    'make_figure_neo2_ql_input_profiles': (
        '.plot_neo2_ql_input_profiles',
        'make_figure_neo2_ql_input_profiles',
    ),
    'add_profiles_to_axes': ('.plot_neo2_ql_input_profiles', 'add_profiles_to_axes'),
    'get_omte_reference_models': ('.plot_omte_reference', 'get_omte_reference_models'),
    'make_figure_omte_reference': ('.plot_omte_reference', 'make_figure_omte_reference'),
    'save_figure_omte_reference': ('.plot_omte_reference', 'save_figure_omte_reference'),
}


def __getattr__(name):
    if name not in _LAZY_EXPORTS:
        raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

    module_name, attr_name = _LAZY_EXPORTS[name]
    module = import_module(module_name, __name__)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value
