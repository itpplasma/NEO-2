"""Plot Om_tE model comparisons against the AUG #30835 reference fixture."""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from .compute_omte import (
    C_CGS,
    compute_poloidal_rotation_neoclassical,
    compute_omte_diamagnetic,
    compute_omte_neoclassical_poloidal,
    compute_omte_toroidal_rotation,
    compute_omte_toroidal_rotation_neo2_convention,
)
from .neo2_output_omte import compute_neo2_omte_from_transport_coefficients
from .neo2_output_omte import decompose_neo2_er_transport_terms
from .neo2_output_omte import compute_neo2_omte_from_k_cof_transport_model


FIXTURE = Path(__file__).resolve().parents[2] / 'test' / 'data' / 'omte_reference_aug30835.npz'
TRANSPORT_FIXTURE = (
    Path(__file__).resolve().parents[2] / 'test' / 'data' / 'neo2_ql_axisymmetric_multispecies_out.h5'
)
LEVEL2_K_SWEEP_MIN = -2.5
LEVEL2_K_SWEEP_MAX = 1.5
LEVEL2_K_SWEEP_COUNT = 161
LEVEL25_DEFAULT_D31_HAT = -1.0
LEVEL25_DEFAULT_K_COF = 0.565


def get_level2_k_scan_reference(
    fixture=FIXTURE,
    k_values=None,
):
    """Return a full-geometry Level 2 k sweep against the stored AUG fixture."""
    ref = np.load(fixture)
    ion_idx = np.where(ref['species_tag'] == ref['species_tag_Vphi'])[0][0]
    z_i = ref['species_def'][0, ion_idx, 0]

    if k_values is None:
        k_values = np.linspace(LEVEL2_K_SWEEP_MIN, LEVEL2_K_SWEEP_MAX, LEVEL2_K_SWEEP_COUNT)
    else:
        k_values = np.asarray(k_values, dtype=float)
    if k_values.ndim != 1 or k_values.size == 0:
        raise ValueError('k_values must be a non-empty 1D array')

    om_scan = np.empty((k_values.size, ref['boozer_s'].size), dtype=float)
    er_scan = np.empty_like(om_scan)
    lvl2_common = {
        'n': ref['n_prof'][:, ion_idx],
        'T': ref['T_prof'][:, ion_idx],
        'dn_ds': ref['dn_ov_ds_prof'][:, ion_idx],
        'dT_ds': ref['dT_ov_ds_prof'][:, ion_idx],
        'z': z_i,
        'aiota': ref['aiota'],
        'sqrtg_bctrvr_phi': ref['sqrtg_bctrvr_phi'],
        'av_nabla_stor': ref['av_nabla_stor'],
        'v_phi': ref['Vphi'],
        'b_theta': ref['bcovar_tht'],
        'b_phi': ref['bcovar_phi'],
    }
    for index, k_value in enumerate(k_values):
        om_scan[index], er_scan[index] = compute_omte_neoclassical_poloidal(
            **lvl2_common,
            k_i=np.full_like(ref['boozer_s'], k_value, dtype=float),
        )

    om_neo2 = C_CGS * ref['Er_neo2'] / (ref['aiota'] * ref['sqrtg_bctrvr_phi'])
    mae = np.mean(np.abs(om_scan - om_neo2), axis=1)
    best_index = int(np.argmin(mae))
    om_min = np.min(om_scan, axis=0)
    om_max = np.max(om_scan, axis=0)
    er_min = np.min(er_scan, axis=0)
    er_max = np.max(er_scan, axis=0)
    q = 1.0 / ref['aiota']
    sqrtg_bctrvr_tht = ref['aiota'] * ref['sqrtg_bctrvr_phi']

    return {
        'boozer_s': ref['boozer_s'],
        'k_values': k_values,
        'om_scan': om_scan,
        'er_scan': er_scan,
        'om_neo2': om_neo2,
        'om_min': om_min,
        'om_max': om_max,
        'er_min': er_min,
        'er_max': er_max,
        'best_k': float(k_values[best_index]),
        'best_mae': float(mae[best_index]),
        'om_best': om_scan[best_index],
        'er_best': er_scan[best_index],
        'neo2_inside_envelope': np.logical_and(om_neo2 >= om_min, om_neo2 <= om_max),
        'q': q,
        'sqrtg_bctrvr_tht': sqrtg_bctrvr_tht,
        'psi_tor_prime_over_reff': np.asarray(ref['sqrtg_bctrvr_phi']),
        'psi_pol_prime_over_reff': sqrtg_bctrvr_tht,
        'psi_tor_prime_over_reff_ov_q': np.asarray(ref['sqrtg_bctrvr_phi']) / q,
    }


def get_omte_reference_models(fixture=FIXTURE):
    """Return Om_tE model curves for the stored AUG reference fixture.

    Level 2 is computed for three collisionality regimes (banana, plateau,
    Pfirsch-Schluter) with K_i = -1.17, -0.5, +0.5 respectively.
    """
    ref = np.load(fixture)
    ion_idx = np.where(ref['species_tag'] == ref['species_tag_Vphi'])[0][0]
    z_i = ref['species_def'][0, ion_idx, 0]

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
    lvl2_common = {
        **common,
        'v_phi': ref['Vphi'],
        'b_theta': ref['bcovar_tht'],
        'b_phi': ref['bcovar_phi'],
    }

    om_neo2 = C_CGS * ref['Er_neo2'] / (ref['aiota'] * ref['sqrtg_bctrvr_phi'])
    om_transport_replay = np.empty_like(ref['boozer_s'], dtype=float)
    er_transport_replay = np.empty_like(ref['boozer_s'], dtype=float)
    for surface_index in range(ref['boozer_s'].size):
        om_transport_replay[surface_index], er_transport_replay[surface_index] = (
            compute_neo2_omte_from_transport_coefficients(
                n_spec=ref['n_prof'][surface_index],
                T_spec=ref['T_prof'][surface_index],
                dn_spec_ov_ds=ref['dn_ov_ds_prof'][surface_index],
                dT_spec_ov_ds=ref['dT_ov_ds_prof'][surface_index],
                species_tag=ref['species_tag'],
                species_tag_vphi=ref['species_tag_Vphi'],
                z_spec=ref['z_spec'],
                Vphi=ref['Vphi'][surface_index],
                aiota=ref['aiota'][surface_index],
                sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'][surface_index],
                av_nabla_stor=ref['av_nabla_stor'][surface_index],
                bcovar_tht=ref['bcovar_tht'][surface_index],
                bcovar_phi=ref['bcovar_phi'][surface_index],
                row_ind=ref['row_ind_spec'][surface_index],
                col_ind=ref['col_ind_spec'][surface_index],
                D31_AX=ref['D31_AX'][surface_index],
                D32_AX=ref['D32_AX'][surface_index],
                D33_AX=ref['D33_AX'][surface_index],
                avEparB_ov_avb2=ref['avEparB_ov_avb2'][surface_index],
                isw_Vphi_loc=0,
            )
        )
    lvl2_k_sweep = get_level2_k_scan_reference(fixture=fixture)
    om_lvl0, er_dia = compute_omte_diamagnetic(**common)
    om_lvl1, er_lvl1 = compute_omte_toroidal_rotation(
        **common,
        v_phi=ref['Vphi'],
        b_theta=ref['bcovar_tht'],
    )

    k_regimes = {'banana': -1.17, 'plateau': -0.5, 'pfirsch_schluter': 0.5}
    om_lvl2 = {}
    er_lvl2 = {}
    for regime, k_val in k_regimes.items():
        k_arr = np.full_like(ref['boozer_s'], k_val, dtype=float)
        om, er = compute_omte_neoclassical_poloidal(**lvl2_common, k_i=k_arr)
        om_lvl2[regime] = om
        er_lvl2[regime] = er

    om_lvl25 = np.empty_like(ref['boozer_s'], dtype=float)
    er_lvl25 = np.empty_like(ref['boozer_s'], dtype=float)
    for surface_index in range(ref['boozer_s'].size):
        om_lvl25[surface_index], er_lvl25[surface_index] = compute_neo2_omte_from_k_cof_transport_model(
            n_spec=ref['n_prof'][surface_index],
            T_spec=ref['T_prof'][surface_index],
            dn_spec_ov_ds=ref['dn_ov_ds_prof'][surface_index],
            dT_spec_ov_ds=ref['dT_ov_ds_prof'][surface_index],
            species_tag=ref['species_tag'],
            species_tag_vphi=ref['species_tag_Vphi'],
            z_spec=ref['z_spec'],
            Vphi=ref['Vphi'][surface_index],
            aiota=ref['aiota'][surface_index],
            sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'][surface_index],
            av_nabla_stor=ref['av_nabla_stor'][surface_index],
            bcovar_tht=ref['bcovar_tht'][surface_index],
            bcovar_phi=ref['bcovar_phi'][surface_index],
            d31_hat=LEVEL25_DEFAULT_D31_HAT,
            k_cof=LEVEL25_DEFAULT_K_COF,
        )

    om_exact, _ = compute_omte_toroidal_rotation_neo2_convention(
        **common,
        vphi=ref['Vphi'],
        bcovar_tht=ref['bcovar_tht'],
        bcovar_phi=ref['bcovar_phi'],
    )
    v_theta = compute_poloidal_rotation_neoclassical(
        dT_ds=ref['dT_ov_ds_prof'][:, ion_idx],
        z=z_i,
        b_phi=ref['bcovar_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
        k_i=np.full_like(ref['boozer_s'], -1.17, dtype=float),
    )
    er_tor = np.asarray(ref['Vphi']) * np.asarray(ref['bcovar_tht']) / C_CGS
    er_pol = -v_theta * np.asarray(ref['bcovar_phi']) / C_CGS

    return {
        'boozer_s': ref['boozer_s'],
        'om_neo2': om_neo2,
        'om_transport_replay': om_transport_replay,
        'om_lvl0': om_lvl0,
        'om_lvl1': om_lvl1,
        'om_lvl2': om_lvl2,
        'lvl2_k_sweep': lvl2_k_sweep,
        'om_lvl25': om_lvl25,
        'om_exact': om_exact,
        'er_neo2': np.asarray(ref['Er_neo2']),
        'er_transport_replay': er_transport_replay,
        'er_dia': er_dia,
        'er_tor': er_tor,
        'er_pol': er_pol,
        'er_lvl1': er_lvl1,
        'er_lvl2': er_lvl2,
        'er_lvl25': er_lvl25,
        'lvl25_parameters': {
            'd31_hat': LEVEL25_DEFAULT_D31_HAT,
            'k_cof': LEVEL25_DEFAULT_K_COF,
        },
    }


def get_transport_reference_decomposition(fixture=TRANSPORT_FIXTURE):
    """Return a full transport-term decomposition for the stored HDF5 fixture."""
    with h5py.File(fixture, 'r') as handle:
        terms = decompose_neo2_er_transport_terms(
            n_spec=np.asarray(handle['n_spec']),
            T_spec=np.asarray(handle['T_spec']),
            dn_spec_ov_ds=np.asarray(handle['dn_spec_ov_ds']),
            dT_spec_ov_ds=np.asarray(handle['dT_spec_ov_ds']),
            species_tag=np.asarray(handle['species_tag']),
            species_tag_vphi=np.asarray(handle['species_tag_Vphi']).reshape(-1)[0],
            z_spec=np.asarray(handle['z_spec']),
            Vphi=float(np.asarray(handle['Vphi'])),
            aiota=float(np.asarray(handle['aiota'])),
            sqrtg_bctrvr_phi=float(np.asarray(handle['sqrtg_bctrvr_phi'])),
            av_nabla_stor=float(np.asarray(handle['av_nabla_stor'])),
            bcovar_tht=float(np.asarray(handle['bcovar_tht'])),
            bcovar_phi=float(np.asarray(handle['bcovar_phi'])),
            row_ind=np.asarray(handle['row_ind_spec']),
            col_ind=np.asarray(handle['col_ind_spec']),
            D31_AX=np.asarray(handle['D31_AX']),
            D32_AX=np.asarray(handle['D32_AX']),
            D33_AX=np.asarray(handle['D33_AX']),
            avEparB_ov_avb2=float(np.asarray(handle['avEparB_ov_avb2'])),
            isw_Vphi_loc=int(np.asarray(handle['isw_Vphi_loc']).reshape(-1)[0]),
        )
        terms['er_stored'] = float(np.asarray(handle['Er']))
    return terms


def _plot_transport_waterfall(axis, terms):
    """Show how the full NEO-2 transport terms move Er away from the common core."""
    labels = ['pressure', 'Vphi', 'D31', 'D32', 'D33', 'denom']
    deltas = np.array(
        [
            terms['er_dia'],
            terms['er_vphi'],
            terms['er_d31'],
            terms['er_d32'],
            terms['er_d33'],
            terms['er_denom'],
        ]
    )
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']
    running = 0.0
    for xpos, (label, delta, color) in enumerate(zip(labels, deltas, colors)):
        axis.bar(xpos, delta, bottom=running, width=0.7, color=color, label=label)
        running += delta
    axis.bar(len(labels), terms['er_total'], width=0.7, color='0.2', label='reconstructed Er')
    axis.axhline(0.0, color='0.7', linewidth=1.0)
    axis.axhline(terms['er_stored'], color='0.1', linewidth=1.2, linestyle='--', label='stored Er')
    axis.set_xticks(range(len(labels) + 1), labels + ['total'])
    axis.set_ylabel('Er [statV/cm]')
    axis.set_title('Full transport waterfall from the exact HDF5 output fixture')
    axis.legend(ncol=4, fontsize=8)


def get_aug_transport_reference_decomposition(fixture=FIXTURE):
    """Return exact AUG transport-term decompositions for all stored surfaces."""
    ref = np.load(fixture)
    terms_by_surface = []
    for surface_index in range(ref['boozer_s'].size):
        terms = decompose_neo2_er_transport_terms(
            n_spec=ref['n_prof'][surface_index],
            T_spec=ref['T_prof'][surface_index],
            dn_spec_ov_ds=ref['dn_ov_ds_prof'][surface_index],
            dT_spec_ov_ds=ref['dT_ov_ds_prof'][surface_index],
            species_tag=ref['species_tag'],
            species_tag_vphi=ref['species_tag_Vphi'],
            z_spec=ref['z_spec'],
            Vphi=ref['Vphi'][surface_index],
            aiota=ref['aiota'][surface_index],
            sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'][surface_index],
            av_nabla_stor=ref['av_nabla_stor'][surface_index],
            bcovar_tht=ref['bcovar_tht'][surface_index],
            bcovar_phi=ref['bcovar_phi'][surface_index],
            row_ind=ref['row_ind_spec'][surface_index],
            col_ind=ref['col_ind_spec'][surface_index],
            D31_AX=ref['D31_AX'][surface_index],
            D32_AX=ref['D32_AX'][surface_index],
            D33_AX=ref['D33_AX'][surface_index],
            avEparB_ov_avb2=ref['avEparB_ov_avb2'][surface_index],
            isw_Vphi_loc=0,
        )
        terms['boozer_s'] = ref['boozer_s'][surface_index]
        terms['er_stored'] = ref['Er_neo2'][surface_index]
        terms_by_surface.append(terms)
    return terms_by_surface


def _plot_aug_transport_terms(axis, terms_by_surface):
    labels = ['pressure', 'Vphi', 'D31', 'D32', 'D33', 'denom', 'total']
    x = np.arange(len(labels))
    width = 0.35
    left = terms_by_surface[0]
    right = terms_by_surface[1]
    values_left = [
        left['er_dia'], left['er_vphi'], left['er_d31'], left['er_d32'],
        left['er_d33'], left['er_denom'], left['er_total'],
    ]
    values_right = [
        right['er_dia'], right['er_vphi'], right['er_d31'], right['er_d32'],
        right['er_d33'], right['er_denom'], right['er_total'],
    ]
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', '0.2']
    for index, (label, color) in enumerate(zip(labels, colors)):
        axis.bar(x[index] - width / 2, values_left[index], width=width, color=color)
        axis.bar(x[index] + width / 2, values_right[index], width=width, color=color)
    axis.axhline(left['er_stored'], color='0.2', linewidth=1.0, linestyle='--')
    axis.axhline(right['er_stored'], color='0.4', linewidth=1.0, linestyle=':')
    axis.axhline(0.0, color='0.7', linewidth=1.0)
    axis.set_xticks(x, labels)
    axis.set_ylabel('Er [statV/cm]')
    axis.set_title('AUG exact transport replay from the rebuilt fixture')
    axis.text(
        0.01,
        0.98,
        f's={left["boozer_s"]:.3f} dashed, s={right["boozer_s"]:.3f} dotted',
        transform=axis.transAxes,
        va='top',
    )


def make_figure_omte_reference(fixture=FIXTURE, transport_fixture=TRANSPORT_FIXTURE):
    """Create a three-panel figure for reduced curves and transport-term breakdown."""
    models = get_omte_reference_models(fixture=fixture)
    transport_terms = get_aug_transport_reference_decomposition(fixture=fixture)
    boozer_s = models['boozer_s']
    lvl2_k_sweep = models['lvl2_k_sweep']

    fig, axes = plt.subplots(3, 1, figsize=(10, 12), constrained_layout=True)

    axes[0].fill_between(
        boozer_s,
        lvl2_k_sweep['om_min'] / 1.0e3,
        lvl2_k_sweep['om_max'] / 1.0e3,
        color='0.88',
        alpha=1.0,
        label='Level 2 k envelope (no geom. corr.)',
        zorder=0,
    )
    axes[0].plot(boozer_s, models['om_neo2'] / 1.0e3, 'o-k', label='NEO-2')
    axes[0].plot(
        boozer_s,
        models['om_transport_replay'] / 1.0e3,
        '-k',
        linewidth=1.0,
        alpha=0.8,
        label='Level 3a replay',
    )
    axes[0].plot(boozer_s, models['om_lvl0'] / 1.0e3, 's--', color='tab:blue', label='Level 0')
    axes[0].plot(
        boozer_s,
        models['om_lvl1'] / 1.0e3,
        'd-.',
        color='tab:orange',
        label='Level 1',
    )
    lvl2_styles = {
        'banana': ('^-', 'tab:green', '$K_i = -1.17$'),
        'plateau': ('v--', 'tab:red', '$K_i = -0.5$'),
        'pfirsch_schluter': ('p:', 'tab:purple', '$K_i = +0.5$'),
    }
    for regime, (fmt, color, label) in lvl2_styles.items():
        axes[0].plot(
            boozer_s,
            models['om_lvl2'][regime] / 1.0e3,
            fmt,
            color=color,
            label=f'Level 2 ({label})',
        )
    axes[0].plot(
        boozer_s,
        lvl2_k_sweep['om_best'] / 1.0e3,
        '--',
        color='0.35',
        linewidth=1.2,
        label=f'Level 2 best k ({lvl2_k_sweep["best_k"]:.2f}, no geom. corr.)',
    )
    lvl25_params = models['lvl25_parameters']
    axes[0].plot(
        boozer_s,
        models['om_lvl25'] / 1.0e3,
        'x-',
        color='tab:brown',
        label=(
            'Level 2.5 '
            f'(D31hat={lvl25_params["d31_hat"]:.3g}, '
            f'k={lvl25_params["k_cof"]:.3g})'
        ),
    )
    axes[0].axhline(0.0, color='0.7', linewidth=1.0)
    axes[0].set_ylabel('Omega_tE [krad/s]')
    axes[0].set_title('Reduced force-balance hierarchy')
    axes[0].legend(fontsize=8)

    x = np.arange(len(boozer_s))
    width = 0.18
    axes[1].bar(x - 1.5 * width, models['er_dia'], width, color='tab:blue', label='pressure')
    axes[1].bar(x - 0.5 * width, models['er_tor'], width, color='tab:orange', label='Vphi')
    axes[1].bar(x + 0.5 * width, models['er_pol'], width, color='tab:green', label='Vtheta')
    axes[1].bar(
        x + 1.5 * width,
        models['er_neo2'],
        width,
        color='0.2',
        label='NEO-2 total',
    )
    axes[1].axhline(0.0, color='0.7', linewidth=1.0)
    axes[1].set_xticks(x, [f's={value:.3f}' for value in boozer_s])
    axes[1].set_ylabel('Er [statV/cm]')
    axes[1].set_title('AUG common terms: what the reduced models actually contain')
    axes[1].legend()

    _plot_aug_transport_terms(axes[2], transport_terms)

    return fig, axes, {'aug': models, 'transport': transport_terms}


def save_figure_omte_reference(
    output_path,
    fixture=FIXTURE,
    transport_fixture=TRANSPORT_FIXTURE,
):
    """Save the Om_tE comparison figure and return the resolved output path."""
    fig, _, _ = make_figure_omte_reference(
        fixture=fixture,
        transport_fixture=transport_fixture,
    )
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
