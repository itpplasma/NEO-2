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
            ion_tag = species_tag_vphi
            for j in range(row.size):
                if int(row.flat[j]) == ion_tag and int(col.flat[j]) == ion_tag:
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


def load_level2_k_band(benchmark_dir=BENCHMARK_DIR,
                       multispec_input=MULTISPEC_INPUT):
    """Compute Level 2 Er for the banana and PS limits and the NEO-2 k_ii."""
    data = load_benchmark_data(benchmark_dir, multispec_input)
    s = data['s']
    er_lvl1 = data['er_lvl1']
    aiota = data['aiota']
    sqrtg_phi = data['sqrtg_bctrvr_phi']

    benchmark_dir = Path(benchmark_dir)
    multispec_input = Path(multispec_input)

    with h5py.File(multispec_input, 'r') as f:
        species_tag = np.asarray(f['species_tag'])
        species_tag_vphi = int(np.asarray(f['species_tag_Vphi']).flat[0])
        ion_idx = int(np.where(species_tag == species_tag_vphi)[0][0])
        z_i = float(np.asarray(f['species_def'])[0, ion_idx, 0])
        boozer_s_in = np.asarray(f['boozer_s'])
        dT_prof = np.asarray(f['dT_ov_ds_prof'])[ion_idx]

    import glob as _glob
    surface_dirs = sorted(_glob.glob(str(benchmark_dir / 'es_*')))
    av_nabla = np.empty(len(surface_dirs))
    for i, d in enumerate(surface_dirs):
        with h5py.File(f'{d}/fulltransp.h5', 'r') as f:
            av_nabla[i] = float(np.asarray(f['avnabpsi']))

    idx_map = np.array([int(np.argmin(np.abs(boozer_s_in - si))) for si in s])
    dT_dr = dT_prof[idx_map] * av_nabla

    k_banana = -1.17
    k_ps = 0.5
    kdg_term_banana = -k_banana * dT_dr / (z_i * E_CGS)
    kdg_term_ps = -k_ps * dT_dr / (z_i * E_CGS)
    kdg_term_kii = -data['k_ii'] * dT_dr / (z_i * E_CGS)

    er_banana = er_lvl1 + kdg_term_banana
    er_ps = er_lvl1 + kdg_term_ps
    er_kii = data['er_lvl2']

    return {
        's': s,
        'er_neo2': data['er_neo2'],
        'er_lvl1': er_lvl1,
        'er_banana': er_banana,
        'er_ps': er_ps,
        'er_kii': er_kii,
        'k_ii': data['k_ii'],
        'aiota': aiota,
        'sqrtg_bctrvr_phi': sqrtg_phi,
    }


def make_figure(data=None, output_path=None):
    """Create the Om_tE figure with Level 1, Level 2 band, and NEO-2."""
    if data is None:
        data = load_level2_k_band()

    s = data['s']
    aiota = data['aiota']
    sqrtg_phi = data['sqrtg_bctrvr_phi']

    def to_om(er):
        return C_CGS * er / (aiota * sqrtg_phi)

    om_neo2 = to_om(data['er_neo2'])
    om_lvl1 = to_om(data['er_lvl1'])
    om_banana = to_om(data['er_banana'])
    om_ps = to_om(data['er_ps'])
    om_kii = to_om(data['er_kii'])

    om_band_lo = np.minimum(om_banana, om_ps)
    om_band_hi = np.maximum(om_banana, om_ps)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    ax.axhline(0, color='gray', lw=0.5)

    ax.fill_between(s, om_band_lo, om_band_hi, alpha=0.2, color='green',
                    label='Level 2 band ($k = -1.17$ to $+0.5$)')

    ax.plot(s, om_banana, 'g-', lw=0.8, alpha=0.5)
    ax.plot(s, om_ps, 'g-', lw=0.8, alpha=0.5)
    ax.plot(s, om_kii, 'g-', lw=2,
            label='Level 2 ($k_{ii}$ from $D_{31}, D_{32}$)')

    ax.plot(s, om_lvl1, 'c-', lw=1.5, label='Level 1 (rotation only)')
    ax.plot(s, om_neo2, 'k-', lw=2.5,
            label='NEO-2 full transport closure', zorder=10)

    ax.set_xlabel('Boozer $s$ (normalised toroidal flux)', fontsize=11)
    ax.set_ylabel(r'$\Omega_{tE}$ [rad/s]', fontsize=11)
    ax.set_xlim(0, 1)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_title('AUG 30835: $\\Omega_{tE}$ from force-balance hierarchy',
                 fontsize=13)
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=200)
    return fig


def load_decomposition_data(benchmark_dir=BENCHMARK_DIR,
                            multispec_input=MULTISPEC_INPUT):
    """Load velocities and Er term decomposition for all surfaces."""
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

    import glob as _glob
    surface_dirs = sorted(_glob.glob(str(benchmark_dir / 'es_*')))
    num = len(surface_dirs)

    s = np.empty(num)
    aiota = np.empty(num)
    sqrtg_phi = np.empty(num)
    bcovar_tht = np.empty(num)
    bcovar_phi = np.empty(num)
    av_nabla = np.empty(num)
    d31_ii = np.empty(num)
    d32_ii = np.empty(num)
    Vtht_neo2 = np.empty(num)
    Vphi_neo2 = np.empty(num)

    for i, d in enumerate(surface_dirs):
        with h5py.File(f'{d}/neo2_multispecies_out.h5', 'r') as f:
            s[i] = float(np.asarray(f['boozer_s']))
            aiota[i] = float(np.asarray(f['aiota']))
            sqrtg_phi[i] = float(np.asarray(f['sqrtg_bctrvr_phi']))
            bcovar_tht[i] = float(np.asarray(f['bcovar_tht']))
            bcovar_phi[i] = float(np.asarray(f['bcovar_phi']))
            D31 = np.asarray(f['D31_AX'])
            D32 = np.asarray(f['D32_AX'])
            row = np.asarray(f['row_ind_spec'])
            col = np.asarray(f['col_ind_spec'])
            ion_tag = species_tag_vphi
            for j in range(row.size):
                if int(row.flat[j]) == ion_tag and int(col.flat[j]) == ion_tag:
                    d31_ii[i] = float(D31.flat[j])
                    d32_ii[i] = float(D32.flat[j])
                    break
            Vtht_neo2[i] = float(np.asarray(f['VthtB_spec']).flat[ion_idx])
            Vphi_neo2[i] = float(np.asarray(f['VphiB_spec']).flat[ion_idx])
        with h5py.File(f'{d}/fulltransp.h5', 'r') as f:
            av_nabla[i] = float(np.asarray(f['avnabpsi']))

    idx_map = np.array([int(np.argmin(np.abs(boozer_s_in - si))) for si in s])
    vphi = vphi_in[idx_map]
    dT_dr = dT_prof[idx_map] * av_nabla
    dp_dr = (T_prof[idx_map] * dn_prof[idx_map]
             + n_prof[idx_map] * dT_prof[idx_map]) * av_nabla
    sqrtg_tht = aiota * sqrtg_phi
    k_ii = 2.5 - d32_ii / d31_ii
    Vtht_kdg = k_ii * C_CGS * dT_dr / (z_i * E_CGS * bcovar_phi)

    er_dia = dp_dr / (n_prof[idx_map] * z_i * E_CGS)
    er_rot = sqrtg_tht * vphi / C_CGS
    er_kdg = -k_ii * dT_dr / (z_i * E_CGS)

    return {
        's': s,
        'Vphi_input': vphi,
        'Vphi_neo2': Vphi_neo2,
        'Vtht_neo2': Vtht_neo2,
        'Vtht_kdg': Vtht_kdg,
        'k_ii': k_ii,
        'er_dia': er_dia,
        'er_rot': er_rot,
        'er_kdg': er_kdg,
    }


def make_decomposition_figure(data=None, model_data=None, output_path=None):
    """Create the 4-panel velocity and Er decomposition figure."""
    if data is None:
        data = load_decomposition_data()
    if model_data is None:
        model_data = load_benchmark_data()

    s = data['s']
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ax = axes[0, 0]
    ax.set_title(r'$V^\varphi$ [rad/s]', fontsize=12)
    ax.plot(s, data['Vphi_neo2'], 'k-', lw=2, label=r'NEO-2 $V^\varphi$ (output)')
    ax.plot(s, data['Vphi_input'], 'r--', lw=1.5, label=r'$V^\varphi$ (input)')
    ax.set_xlabel('Boozer $s$')
    ax.set_ylabel(r'$V^\varphi$ [rad/s]')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

    ax = axes[0, 1]
    ax.set_title(r'$V^\vartheta$ [rad/s]', fontsize=12)
    ax.plot(s, data['Vtht_neo2'], 'k-', lw=2, label=r'NEO-2 $V^\vartheta$ (full)')
    ax.plot(s, data['Vtht_kdg'], 'g--', lw=1.5, label=r'KDG estimate ($k_{ii}$)')
    ax.axhline(0, color='gray', lw=0.5)
    ax.set_xlabel('Boozer $s$')
    ax.set_ylabel(r'$V^\vartheta$ [rad/s]')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

    ax = axes[1, 0]
    ax.set_title('$E_r$ term decomposition [statV/cm]', fontsize=12)
    ax.axhline(0, color='gray', lw=0.5)
    ax.fill_between(s, 0, data['er_dia'], alpha=0.3, color='blue',
                    label=r'$\frac{1}{Z_i e n_i}\frac{dp_i}{dr}$')
    ax.fill_between(s, 0, data['er_rot'], alpha=0.3, color='cyan',
                    label=r'$\frac{\sqrt{g}B^\vartheta}{c} V^\varphi$')
    ax.fill_between(s, 0, data['er_kdg'], alpha=0.3, color='green',
                    label=r'$-\frac{k_i}{Z_i e}\frac{dT_i}{dr}$')
    ax.plot(s, data['er_dia'] + data['er_rot'] + data['er_kdg'],
            'g-', lw=2, label='Level 2 total')
    ax.plot(s, model_data['er_neo2'], 'k-', lw=2, label='NEO-2 full')
    ax.set_xlabel('Boozer $s$')
    ax.set_ylabel('$E_r$ contribution [statV/cm]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

    ax = axes[1, 1]
    ax.set_title('$E_r$ model comparison [statV/cm]', fontsize=12)
    ax.axhline(0, color='gray', lw=0.5)
    ax.plot(s, model_data['er_neo2'], 'k-', lw=2, label='NEO-2 full', zorder=10)
    ax.plot(s, model_data['er_lvl2'], 'g-', lw=1.5, label='Level 2 ($k_{ii}$)')
    ax.plot(s, model_data['er_reduced'], 'm--', lw=1.5, label='Reduced single-ion')
    ax.plot(s, model_data['er_lvl1'], 'c-', lw=1, alpha=0.7, label='Level 1')
    ax.plot(s, model_data['er_lvl0'], 'b-', lw=1, alpha=0.7, label='Level 0')
    ax.set_xlabel('Boozer $s$')
    ax.set_ylabel('$E_r$ [statV/cm]')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

    fig.suptitle('AUG 30835: Velocities, $E_r$ decomposition, and model hierarchy',
                 fontsize=14, fontweight='bold')
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=200)
    return fig


if __name__ == '__main__':
    band_data = load_level2_k_band()
    make_figure(band_data, '/tmp/aug30835_omte.png')
    print('Saved /tmp/aug30835_omte.png')

    model_data = load_benchmark_data()
    decomp_data = load_decomposition_data()
    make_decomposition_figure(decomp_data, model_data,
                              '/tmp/aug30835_decomposition.png')
    print('Saved /tmp/aug30835_decomposition.png')
