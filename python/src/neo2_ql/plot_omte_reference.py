"""Plot Om_tE model comparisons against the AUG #30835 reference fixture."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .compute_omte import (
    C_CGS,
    compute_omte_diamagnetic,
    compute_omte_neoclassical_poloidal_auto_k,
    compute_omte_toroidal_rotation,
    compute_omte_toroidal_rotation_neo2_convention,
)


FIXTURE = Path(__file__).resolve().parents[2] / 'test' / 'data' / 'omte_reference_aug30835.npz'


def get_omte_reference_models(fixture=FIXTURE, nu_star=None):
    """Return Om_tE model curves for the stored AUG reference fixture."""
    ref = np.load(fixture)
    ion_idx = np.where(ref['species_tag'] == ref['species_tag_Vphi'])[0][0]
    z_i = ref['species_def'][0, ion_idx, 0]

    if nu_star is None:
        nu_star = np.full_like(ref['boozer_s'], 0.05, dtype=float)
    else:
        nu_star = np.asarray(nu_star, dtype=float)

    common = {
        'n': ref['n_prof'][:, ion_idx],
        'T': ref['T_prof'][:, ion_idx],
        'dn_ds': ref['dn_ov_ds_prof'][:, ion_idx],
        'dT_ds': ref['dT_ov_ds_prof'][:, ion_idx],
        'z': z_i,
        'aiota': ref['aiota'],
        'sqrtg_bctrvr_phi': ref['sqrtg_bctrvr_phi'],
        'av_nabla_stor': ref['av_nabla_stor'],
    }

    om_neo2 = C_CGS * ref['Er_neo2'] / (ref['aiota'] * ref['sqrtg_bctrvr_phi'])
    om_lvl0, _ = compute_omte_diamagnetic(**common)
    om_lvl1, _ = compute_omte_toroidal_rotation(
        **common,
        v_phi=ref['R0'] * ref['Vphi'],
        b_theta=ref['bcovar_tht'] / ref['R0'],
    )
    om_lvl2, _ = compute_omte_neoclassical_poloidal_auto_k(
        **common,
        v_phi=ref['R0'] * ref['Vphi'],
        b_theta=ref['bcovar_tht'] / ref['R0'],
        b_phi=ref['bcovar_phi'] / ref['R0'],
        nu_star=nu_star,
    )
    om_exact, _ = compute_omte_toroidal_rotation_neo2_convention(
        **common,
        vphi=ref['Vphi'],
        bcovar_tht=ref['bcovar_tht'],
        bcovar_phi=ref['bcovar_phi'],
    )

    return {
        'boozer_s': ref['boozer_s'],
        'om_neo2': om_neo2,
        'om_lvl0': om_lvl0,
        'om_lvl1': om_lvl1,
        'om_lvl2': om_lvl2,
        'om_exact': om_exact,
        'nu_star': nu_star,
    }


def make_figure_omte_reference(fixture=FIXTURE, nu_star=None):
    """Create a two-panel figure for practical and strict Om_tE models."""
    models = get_omte_reference_models(fixture=fixture, nu_star=nu_star)
    boozer_s = models['boozer_s']

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True, constrained_layout=True)

    axes[0].plot(boozer_s, models['om_neo2'] / 1.0e3, 'o-k', label='NEO-2')
    axes[0].plot(boozer_s, models['om_lvl0'] / 1.0e3, 's--', color='tab:blue', label='Level 0')
    axes[0].plot(boozer_s, models['om_lvl1'] / 1.0e3, 'd-.', color='tab:orange', label='Level 1 proxy')
    axes[0].plot(
        boozer_s,
        models['om_lvl2'] / 1.0e3,
        '^-',
        color='tab:green',
        label='Level 2 proxy (auto K_i)',
    )
    axes[0].axhline(0.0, color='0.7', linewidth=1.0)
    axes[0].set_ylabel('Omega_tE [krad/s]')
    axes[0].set_title('Practical force-balance models')
    axes[0].legend()

    axes[1].plot(boozer_s, models['om_neo2'] / 1.0e6, 'o-k', label='NEO-2')
    axes[1].plot(
        boozer_s,
        models['om_exact'] / 1.0e6,
        'x-',
        color='tab:red',
        linewidth=2.0,
        label='Strict NEO-2 Vphi convention',
    )
    axes[1].axhline(0.0, color='0.7', linewidth=1.0)
    axes[1].set_xlabel('s_tor')
    axes[1].set_ylabel('Omega_tE [Mrad/s]')
    axes[1].set_title('Reduced isw_Vphi_loc=0 algebra')
    axes[1].legend()

    return fig, axes, models


def save_figure_omte_reference(output_path, fixture=FIXTURE, nu_star=None):
    """Save the Om_tE comparison figure and return the resolved output path."""
    fig, _, _ = make_figure_omte_reference(fixture=fixture, nu_star=nu_star)
    output_path = Path(output_path).resolve()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return output_path


def main():
    """Write the default Om_tE comparison figure to a file."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output',
        default='/tmp/omte_reference_models.png',
        help='Path for the output PNG',
    )
    args = parser.parse_args()
    output_path = save_figure_omte_reference(args.output)
    print(output_path)


if __name__ == '__main__':
    main()
