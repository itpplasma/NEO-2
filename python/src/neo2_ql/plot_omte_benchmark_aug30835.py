"""Plot Om_tE model hierarchy using the full 100-surface AUG 30835 benchmark."""

from pathlib import Path

import glob
import h5py
import matplotlib.pyplot as plt
import numpy as np

from .compute_omte import C_CGS, E_CGS


BENCHMARK_DIR = Path(
    '/home/ert/data/AUG/NEO-2/30835/2016_controlled_fusion_rmp90_benchmark'
)
MULTISPEC_INPUT = BENCHMARK_DIR / 'aug_30835.in'


def load_benchmark_data(benchmark_dir=BENCHMARK_DIR,
                        multispec_input=MULTISPEC_INPUT):
    """Load input profiles and NEO-2 output from all surfaces."""
    benchmark_dir = Path(benchmark_dir)
    multispec_input = Path(multispec_input)

    with h5py.File(multispec_input, 'r') as f:
        boozer_s_in = np.asarray(f['boozer_s'])
        species_tag = np.asarray(f['species_tag'])
        species_tag_vphi = int(np.asarray(f['species_tag_Vphi']).flat[0])
        ion_idx = int(np.where(species_tag == species_tag_vphi)[0][0])
        z_i = float(np.asarray(f['species_def'])[0, ion_idx, 0])
        n_prof = np.asarray(f['n_prof'])[ion_idx]
        T_prof = np.asarray(f['T_prof'])[ion_idx]
        dn_prof = np.asarray(f['dn_ov_ds_prof'])[ion_idx]
        dT_prof = np.asarray(f['dT_ov_ds_prof'])[ion_idx]
        vphi_in = np.asarray(f['Vphi'])

    surface_dirs = sorted(glob.glob(str(benchmark_dir / 'es_*')))
    num_surfaces = len(surface_dirs)

    s = np.empty(num_surfaces)
    er_neo2 = np.empty(num_surfaces)
    aiota = np.empty(num_surfaces)
    sqrtg_phi = np.empty(num_surfaces)
    bcovar_tht = np.empty(num_surfaces)
    bcovar_phi = np.empty(num_surfaces)
    av_nabla = np.empty(num_surfaces)
    d31_ii = np.empty(num_surfaces)
    d32_ii = np.empty(num_surfaces)

    for i, d in enumerate(surface_dirs):
        out_h5 = Path(d) / 'neo2_multispecies_out.h5'
        ft_h5 = Path(d) / 'fulltransp.h5'
        with h5py.File(out_h5, 'r') as f:
            s[i] = float(np.asarray(f['boozer_s']))
            er_neo2[i] = float(np.asarray(f['Er']))
            aiota[i] = float(np.asarray(f['aiota']))
            sqrtg_phi[i] = float(np.asarray(f['sqrtg_bctrvr_phi']))
            bcovar_tht[i] = float(np.asarray(f['bcovar_tht']))
            bcovar_phi[i] = float(np.asarray(f['bcovar_phi']))
            D31 = np.asarray(f['D31_AX'])
            D32 = np.asarray(f['D32_AX'])
            row = np.asarray(f['row_ind_spec'])
            col = np.asarray(f['col_ind_spec'])
            for j in range(row.size):
                if int(row.flat[j]) == ion_idx and int(col.flat[j]) == ion_idx:
                    d31_ii[i] = float(D31.flat[j])
                    d32_ii[i] = float(D32.flat[j])
                    break
        with h5py.File(ft_h5, 'r') as f:
            av_nabla[i] = float(np.asarray(f['avnabpsi']))

    idx_map = np.array([int(np.argmin(np.abs(boozer_s_in - si))) for si in s])
    n = n_prof[idx_map]
    T = T_prof[idx_map]
    dn = dn_prof[idx_map]
    dT = dT_prof[idx_map]
    vphi = vphi_in[idx_map]

    sqrtg_tht = aiota * sqrtg_phi
    k_ii = 2.5 - d32_ii / d31_ii
    dp_dr = (T * dn + n * dT) * av_nabla
    dT_dr = dT * av_nabla

    er_lvl0 = dp_dr / (n * z_i * E_CGS)
    er_lvl1 = er_lvl0 + sqrtg_tht * vphi / C_CGS
    er_lvl2 = er_lvl1 - k_ii * dT_dr / (z_i * E_CGS)
    geom = bcovar_phi / (z_i * E_CGS * (aiota * bcovar_tht + bcovar_phi))
    er_red = sqrtg_tht * vphi / C_CGS + er_lvl0 - k_ii * geom * dT_dr

    return {
        's': s,
        'er_neo2': er_neo2,
        'er_lvl0': er_lvl0,
        'er_lvl1': er_lvl1,
        'er_lvl2': er_lvl2,
        'er_reduced': er_red,
        'k_ii': k_ii,
        'aiota': aiota,
        'sqrtg_bctrvr_phi': sqrtg_phi,
    }


def make_figure(data=None, output_path=None):
    """Create the 2-panel Er / Om_tE comparison figure."""
    if data is None:
        data = load_benchmark_data()

    s = data['s']
    aiota = data['aiota']
    sqrtg_phi = data['sqrtg_bctrvr_phi']

    om = {k: C_CGS * data[k] / (aiota * sqrtg_phi)
          for k in ['er_neo2', 'er_lvl0', 'er_lvl1', 'er_lvl2', 'er_reduced']}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    ax1.set_title('$E_r$ [statV/cm]', fontsize=13)
    ax1.axhline(0, color='gray', lw=0.5)
    ax1.plot(s, data['er_neo2'], 'k-', lw=2, label='NEO-2 full', zorder=10)
    ax1.plot(s, data['er_lvl2'], 'g-', lw=1.5, label='Level 2 ($k_{ii}$)')
    ax1.plot(s, data['er_reduced'], 'm--', lw=1.5,
             label='Reduced single-ion ($k_{ii}$)')
    ax1.plot(s, data['er_lvl1'], 'c-', lw=1, alpha=0.7, label='Level 1')
    ax1.plot(s, data['er_lvl0'], 'b-', lw=1, alpha=0.7, label='Level 0')
    ax1.set_xlabel('Boozer $s$')
    ax1.set_ylabel('$E_r$ [statV/cm]')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)

    ax2.set_title(r'$\Omega_{tE}$ [rad/s]', fontsize=13)
    ax2.axhline(0, color='gray', lw=0.5)
    ax2.plot(s, om['er_neo2'], 'k-', lw=2, label='NEO-2 full', zorder=10)
    ax2.plot(s, om['er_lvl2'], 'g-', lw=1.5, label='Level 2 ($k_{ii}$)')
    ax2.plot(s, om['er_reduced'], 'm--', lw=1.5,
             label='Reduced single-ion ($k_{ii}$)')
    ax2.plot(s, om['er_lvl1'], 'c-', lw=1, alpha=0.7, label='Level 1')
    ax2.plot(s, om['er_lvl0'], 'b-', lw=1, alpha=0.7, label='Level 0')
    ax2.set_xlabel('Boozer $s$')
    ax2.set_ylabel(r'$\Omega_{tE}$ [rad/s]')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 1)

    fig.suptitle('AUG 30835: Force-balance model hierarchy (100 surfaces)',
                 fontsize=13, fontweight='bold')
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=200)
    return fig


if __name__ == '__main__':
    data = load_benchmark_data()
    make_figure(data, '/tmp/aug30835_100surf.png')
    print('Saved /tmp/aug30835_100surf.png')
