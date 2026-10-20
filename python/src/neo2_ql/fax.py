"""Lossless reader and independent oracles for NEO-2-QL f_AX exports.

Fortran indices in offset/sample_row are one based offsets as documented in
DOC/fax_export.md. Arrays are returned in Fortran axis order, not h5py order.
"""
from pathlib import Path
import json
import math
import numpy as np
import h5py


def read_fax(path):
    with h5py.File(path) as f:
        d = {k: np.asarray(v[()]).T.copy() for k, v in f.items()}
    for k in ("format_version", "field_species", "num_species", "lag",
              "collop_base_exp", "isw_lorentz", "isw_calc_Er"):
        d[k] = int(d[k].item())
    if d["format_version"] != 1:
        raise ValueError(f"{path}: unsupported f_AX format")
    n = d["num_species"]
    if d["g_real"].shape[1:] != (3, n):
        raise ValueError(f"{path}: missing source/species axes")
    if abs(float(d['phi_direction'].item())) != 1:
        raise ValueError(f'{path}: invalid parallel-branch orientation')
    if not np.isfinite(d["g_real"]).all():
        raise ValueError(f"{path}: nonfinite distribution")
    if np.max(np.abs(d["g_imag"])) > 1e-12*max(1., np.max(np.abs(d["g_real"]))):
        raise ValueError(f"{path}: axisymmetric solution is complex")
    return d


def basis(d, x):
    x = float(x)
    if x < 0:
        raise ValueError("x=v/v_T must be nonnegative")
    if not d["basis_extend"].item() and x > d["basis_edges"][-1]:
        return np.zeros(d["lag"]+1)
    j = np.clip(np.searchsorted(d["basis_edges"], x, side="right")-1,
                0, len(d["basis_centers"])-1)
    return np.polynomial.polynomial.polyval(x-d["basis_centers"][j],
                                          d["basis_coeff"][:, :, j])


def block(d, step):
    """Band-integrated g: (energy block, signed band, source, drive species)."""
    n = int(d["npassing"][step])+1
    off = int(d["offset"][step])
    nb = d["lag"]+1
    return d["g_real"][off:off+2*n*nb].reshape(nb, 2*n, 3, d["num_species"])


def oracles(d):
    """Native phi_exp samples and independent energy quadrature -> qflux.

    Norms use max(1, |reference|) for distribution samples; transport uses a
    matrix infinity scale per source species, preserving near-zero entries.
    The monoenergetic Lorentz model has no Maxwellian energy integral.
    """
    from scipy.integrate import quad
    basis_error = max(np.max(np.abs(basis(d, x)-d["sample_basis"][:, j]) /
                             np.maximum(1., np.abs(d["sample_basis"][:, j])))
                      for j, x in enumerate(d["sample_x"]))
    point_error = 0.
    for t, st in enumerate(d["sample_step"]):
        n = int(d["npassing"][st])+1
        row = int(d["sample_row"][t])-1
        ids = row+2*n*np.arange(d["lag"]+1)
        for j, x in enumerate(d["sample_x"]):
            got = np.einsum('mkb,m->kb', d["g_real"][ids], basis(d, x))
            ref = d["sample_g"][:, :, j, t]
            point_error = max(point_error, float(np.max(np.abs(got-ref) /
                                                       np.maximum(1., np.abs(ref)))))
    out = dict(basis_error=float(basis_error), pointwise_error=point_error,
               pointwise_tolerance=2e-10, samples=8*6*3*d["num_species"])
    if max(basis_error, point_error) > out["pointwise_tolerance"]:
        raise ValueError(f"native distribution reconstruction failed: {out}")
    if d["isw_lorentz"] == 1:
        out["energy_integral"] = "not defined for monoenergetic Lorentz model"
        return out
    weights = np.empty((3, d["lag"]+1))
    for k, power in enumerate((4, 3, 6)):  # NEO-2 swaps the last two sources
        for m in range(d["lag"]+1):
            edges = list(d["basis_edges"])
            if d["basis_extend"].item():
                edges[-1] = np.inf
            weights[k, m] = sum(quad(lambda x: math.exp(-x*x)*x**power*basis(d,x)[m],
                                      lo, hi, epsabs=1e-11, epsrel=1e-11)[0]
                                  for lo, hi in zip(edges[:-1], edges[1:]))/math.sqrt(math.pi)
    # Build the phase-space convolution using independently integrated phi_m,
    # not the exported NEO-2 weightlag or flux_vector.
    flux = np.empty_like(d["flux_geometry"])
    for st in range(len(d["npassing"])):
        n = int(d["npassing"][st])+1
        off = int(d["offset"][st])
        for m in range(d["lag"]+1):
            sl = slice(off+2*n*m, off+2*n*(m+1))
            flux[:, sl] = d["flux_geometry"][:, sl]*weights[:, m, None]
    got = np.einsum('ji,ikb->jkb', flux, d["g_real"])
    ref = d["qflux"]
    scale = np.maximum(1e-30, np.max(np.abs(ref), axis=(0, 1)))
    err = float(np.max(np.abs(got-ref)/scale))
    entry_scale = np.maximum(np.abs(ref), 1e-12*scale[None,None,:])
    entry_error = float(np.max(np.abs(got-ref)/entry_scale))
    out.update(energy_weight_error=float(np.max(np.abs(weights-d["weightlag"]))),
               transport_error=err, transport_entrywise_error=entry_error,
               transport_tolerance=2e-8,
               entrywise_zero_floor_fraction=1e-12,
               qflux_native=ref.tolist(), qflux_reconstructed=got.tolist())
    if max(err, entry_error) > out["transport_tolerance"]:
        raise ValueError(f"energy-integrated qflux reconstruction failed: {out}")
    return out


def source_amplitudes(d, er, eparb_over_b2=0.):
    """beta_b A_k^b in CODE source order (density, parallel E, temperature).

    DOC/ripple_solver_normalizations_and_output.tex, multispec:fa1norm and
    multispec:fa1basisfun. Er is in statV/cm; all metadata are CGS.
    """
    n, t, m, z = (d[k] for k in ('density_cm3', 'temperature_erg', 'mass_g',
                                 'charge_number'))
    grad = d['y_ntv'][6]/d['y_ntv'][5]
    avbhat = d['y_ntv'][13]/d['y_ntv'][12]
    beta = avbhat*np.sqrt(2*t/m)*m*2.9979e10/(z*4.8032e-10*d['B_ref_G'])/grad
    a2 = grad*d['dT_ds_erg']/t
    a1 = grad*d['dn_ds_cm3']/n-1.5*a2-er*z*4.8032e-10/t
    a3 = z*4.8032e-10/t*eparb_over_b2
    return np.asarray([beta*a1, d['B_ref_G']/avbhat*a3, beta*a2])


def dimensional_oracle(d, native_output, reconstructed_qflux):
    """Independent CGS coefficient check against ntv_mod's public Dij output.

    Keep the native FSA factors beta_out in ntv_mod:286. This additionally
    checks charge signs, field/driving species order and the 1,3,2 source map.
    """
    m,t,z = (d[k] for k in ('mass_g','temperature_erg','charge_number'))
    y = d['y_ntv']
    grad, avbhat = y[6]/y[5], y[13]/y[12]
    rho = np.sqrt(2*t/m)*m*2.9979e10/(z*4.8032e-10*d['B_ref_G'])
    beta = np.stack([avbhat*rho/grad, avbhat*rho/grad,
                     np.full_like(rho,float(d['B_ref_G'])/avbhat)])
    a,n = d['field_species'],d['num_species']
    q = np.asarray(reconstructed_qflux)
    mapping = (0,2,1)
    got, ref = np.empty((3,3,n)), np.empty((3,3,n))
    with h5py.File(native_output) as h:
        for j in range(3):
            for k in range(3):
                got[j,k] = -np.sqrt(2*t[a]/m[a])*beta[j,a]*beta[k]/y[5] * q[mapping[j],mapping[k]]
                ref[j,k] = np.asarray(h[f'D{j+1}{k+1}_AX'][()]).reshape(n,n)[a]
    scale = np.maximum(np.abs(ref),1e-12*np.max(np.abs(ref)))
    error = float(np.max(np.abs(got-ref)/np.maximum(1e-30,scale)))
    result = dict(max_entrywise_error=error,tolerance=2e-8,avbhat=float(avbhat),
                  D11_native_cm2_s=ref[0,0].tolist(),
                  D11_reconstructed_cm2_s=got[0,0].tolist())
    if error > result['tolerance']:
        raise ValueError(f'native dimensional transport normalization failed: {result}')
    return result


def pitch_rows(d, step):
    """Band means of signed lambda and normalized g (m, lambda, k, b)."""
    n = int(d['npassing'][step])+1
    b = d['bhat'][step]
    lo = d['eta'][:n]
    hi = np.r_[d['eta'][1:n], 1./b]
    width = hi-lo
    if np.any(width <= 0):
        raise ValueError(f"step {step}: nonpositive pitch band")
    # Exact band mean of |lambda|, identical to the legacy circular compiler.
    lam = 2./(3.*b)*(np.maximum(1-b*lo,0)**1.5-
                     np.maximum(1-b*hi,0)**1.5)/width
    g = block(d, step)
    # Counter-passing occupies reversed bands; positive lambda increases as eta falls.
    xx = np.r_[-lam, lam[::-1]]
    yy = np.concatenate([g[:, n:][:, ::-1]/width[None,:,None,None],
                         g[:, :n][:, ::-1]/width[None,::-1,None,None]], axis=1)
    if d['phi_direction'].item() < 0:
        xx, yy = -xx[::-1], yy[:, ::-1]
    return xx, yy


def position(d):
    b = d['bhat']
    signs = np.sign(d['dlogb_dphi']).copy()
    nonzero = np.flatnonzero(signs)
    if not len(nonzero):
        raise ValueError('field has no position branch')
    # An extremum has no derivative sign. Use a neighbouring branch there;
    # the two maximum labels +/-pi are the same periodic position.
    for j in np.flatnonzero(signs == 0):
        previous = nonzero[nonzero < j]
        signs[j] = signs[previous[-1] if len(previous) else nonzero[0]]
    return d['phi_direction'].item()*signs*np.arccos(np.clip(
        1-2*(b-b.min())/(b.max()-b.min()),-1,1))


def admission(d):
    """Reject geometry that the single-well vartheta chart cannot represent."""
    b = d['bhat']
    derivative = d['dlogb_dphi']
    sig = np.sign(derivative[np.abs(derivative) > 1e-8*np.max(np.abs(derivative))])
    changes = int(np.count_nonzero(sig[1:] != sig[:-1]))
    if b.max() <= b.min() or len(sig) < 2 or changes > 2:
        raise ValueError(f"field species {d['field_species']}: single-well admission failed; "
                         f"{changes} interior branch changes")
    return dict(interior_branch_changes=changes, Bmin_hat=float(b.min()),
                Bmax_hat=float(b.max()), position_chart='single well, sign(b^phi dlogB/dphi)',
                phi_direction=float(d['phi_direction'].item()),
                energy_basis=int(d['collop_base_exp']), energy_blocks=d['lag']+1,
                driving_species=d['num_species'], field_species=d['field_species'],
                physical_transport_qualified=False)


def regular_blocks(d, theta, lambdas):
    # C order output axes: basis, source, driving species, theta, lambda.
    vt = position(d)
    order = np.argsort(vt, kind='stable')
    keep = np.r_[True, np.diff(vt[order]) > 1e-12]
    order = order[keep]
    v = vt[order]
    v = np.r_[v-2*np.pi, v, v+2*np.pi]
    shape = (d['lag']+1, 3, d['num_species'], len(theta), len(lambdas))
    values, slopes = np.empty(shape), np.empty(shape)
    for m in range(d['lag']+1):
        rows, derivs = [], []
        for st in order:
            x, y = pitch_rows(d, st)
            rows.append(np.stack([[np.interp(lambdas, x, y[m,:,k,b])
                                   for b in range(d['num_species'])] for k in range(3)]))
            dy = np.gradient(y[m], x, axis=0)
            derivs.append(np.stack([[np.interp(lambdas, x, dy[:,k,b])
                                     for b in range(d['num_species'])] for k in range(3)]))
        for out, data in ((values, rows), (slopes, derivs)):
            data = np.asarray(data)
            data = np.concatenate([data, data, data], axis=0)
            for k in range(3):
                for b in range(d['num_species']):
                    for l in range(len(lambdas)):
                        out[m,k,b,:,l] = np.interp(theta, v, data[:,k,b,l])
    # Native consumer order (lambda, theta, basis, source, driving species).
    return values.transpose(4,3,0,1,2), slopes.transpose(4,3,0,1,2)


def regular_transport_oracle(d, values, theta, lambdas):
    """Energy-integrate the actual regular table, then use native band quadrature.

    Linearity lets the independently verified energy weights integrate each
    polynomial block before position/pitch interpolation. The source forces
    are kept separate: this checks D11, not a cancellation of profile forces.
    """
    integrated = np.einsum('ltmkb,m->ltkb', values, d['weightlag'][0])
    got = np.zeros((3,d['num_species']))
    vt = position(d)
    for st in range(len(d['npassing'])):
        n = int(d['npassing'][st])+1
        x, _ = pitch_rows(d,st)
        it = int(np.clip(np.searchsorted(theta,vt[st])-1,0,len(theta)-2))
        ft = np.clip((vt[st]-theta[it])/(theta[it+1]-theta[it]),0,1)
        row = (1-ft)*integrated[:,it]+ft*integrated[:,it+1]
        il = np.clip(np.searchsorted(lambdas,x)-1,0,len(lambdas)-2)
        fl = np.clip((x-lambdas[il])/(lambdas[il+1]-lambdas[il]),0,1)
        band = (1-fl)[:,None,None]*row[il]+fl[:,None,None]*row[il+1]
        if d['phi_direction'].item() < 0:
            band = band[::-1]
        native = np.concatenate([band[n:][::-1],band[:n][::-1]])
        hi = np.r_[d['eta'][1:n],1/d['bhat'][st]]
        widths = hi-d['eta'][:n]
        native *= np.r_[widths,widths[::-1]][:,None,None]
        off = int(d['offset'][st])
        got += np.einsum('i,ikb->kb',d['flux_geometry'][0,off:off+2*n],native)
    # Driving-species rho factors convert the complete density-drive row to
    # physical D11; the common field-species factor cancels in this ratio.
    rho = np.sqrt(2*d['temperature_erg']/d['mass_g'])*d['mass_g']/d['charge_number']
    ref = d['qflux'][0,0]*rho
    value = got[0]*rho
    error = float(np.max(np.abs(value-ref))/max(1e-30,np.max(np.abs(ref))))
    result = dict(D11_row_error=error,tolerance=0.02,
                  norm='native D11 field-row maximum, all driving species',
                  qflux_density_native=d['qflux'][0,0].tolist(),
                  qflux_density_regular_table=got[0].tolist())
    if error > result['tolerance']:
        raise ValueError(f'regular-table energy-integrated D11 failed: {result}')
    return result


def write_table(exports, path, er, epar=0., nt=181, nl=401):
    import hashlib
    ds = sorted((read_fax(f) for f in exports), key=lambda d: d['field_species'])
    d = ds[0]
    n = d['num_species']
    if len(ds) != n or [v['field_species'] for v in ds] != list(range(n)):
        raise ValueError('one export for EACH field species is required')
    for other in ds[1:]:
        for k in ('species_tag', 'mass_g', 'charge_number', 'temperature_erg',
                  'density_cm3', 'basis_coeff', 'basis_edges', 'basis_extend', 's_tor'):
            if not np.array_equal(d[k], other[k]):
                raise ValueError(f'field species metadata differ: {k}')
    checks = [admission(v) for v in ds]
    amp = source_amplitudes(d, er, epar)
    theta, lambdas = np.linspace(-np.pi,np.pi,nt), np.linspace(-1,1,nl)
    degree, nb, segments = d['basis_coeff'].shape
    records = [d['mass_g'], d['charge_number'], d['temperature_erg'], amp,
               d['basis_edges'], d['basis_centers'], d['basis_coeff'], theta, lambdas]
    blocks = [regular_blocks(v,theta,lambdas) for v in ds]
    transport = [regular_transport_oracle(v,b[0],theta,lambdas) for v,b in zip(ds,blocks)]
    # Retain each source independently; no profile-gradient division or lost cross-species term.
    records += [np.stack([b[channel] for b in blocks],axis=-1) for channel in (0,1)]
    path = Path(path)
    with path.open('w') as f:
        f.write('# NEO-2 physical f_AX/f_M; source order density, E_parallel, temperature\n')
        f.write('FAX_ENERGY_V1\n')
        f.write(f'{nt} {nl} {nb} {n} {degree-1} {segments} '
                f'{"T" if d["basis_extend"].item() else "F"}\n')
        for r in records:
            flat = np.asarray(r).ravel(order='F')
            # Bounded records also keep Fortran list-directed input portable.
            for off in range(0,len(flat),6):
                f.write(' '.join(f'{float(x):.17e}' for x in flat[off:off+6])+'\n')
    receipt = dict(regular_table_transport=transport, table=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                   source_exports={str(f): hashlib.sha256(Path(f).read_bytes()).hexdigest()
                                   for f in exports}, source_amplitudes=amp.tolist(),
                   Er_statV_per_cm=er, EparB_over_B2=epar, admission=checks,
                   dimensions=dict(theta=nt,pitch=nl,basis=nb,species=n),
                   transport_qualified=False)
    path.with_suffix(path.suffix+'.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def consumer_oracle(d, table, executable, amplitudes):
    """Compare native distribution band moments with the actual APPLETS reader."""
    import subprocess
    if not np.any(amplitudes):
        raise ValueError('consumer oracle requires nonzero source amplitudes')
    maxerr = 0.
    samples, queries, refs = [], [], []
    pitches = []
    for sample in range(8):
        st = int(d['sample_step'][sample])
        n = int(d['npassing'][st])+1
        x, values = pitch_rows(d,st)
        tp = np.searchsorted(d['eta'],1/np.max(d['bhat']),side='right')-1
        # Passing, trapped-passing boundary and the incomplete turning band.
        bands = sorted({0, int(d['sample_band'][sample])-1, int(np.clip(tp,0,n-1)), n-1})
        for band in bands:
            col = band if d['sample_sign'][sample] < 0 else 2*n-band-1
            if d['phi_direction'].item() < 0:
                col = 2*n-col-1
            theta, pitch = float(position(d)[st]), float(x[col])
            pitches.append(pitch)
            lo = d['eta'][band]
            hi = d['eta'][band+1] if band+1 < n else 1/d['bhat'][st]
            row = int(d['offset'][st])+band
            if d['sample_sign'][sample] < 0:
                row = int(d['offset'][st])+2*n-band-1
            ids = row+2*n*np.arange(d['lag']+1)
            for j in (1,2,3,4,5):
                energy = float(d['temperature_erg'][d['field_species']]*d['sample_x'][j]**2)
                queries.append((d['field_species']+1, theta, pitch, energy))
                # Native phi_exp values (not the exported polynomial evaluator)
                # and native band integrals provide the independent moment.
                ref = float(np.einsum('mkb,m,kb->',d['g_real'][ids],
                                     d['sample_basis'][:,j],amplitudes)/(hi-lo))
                refs.append(ref)
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat') as f:
        np.savetxt(f, queries, fmt=['%d','%.17e','%.17e','%.17e'])
        f.flush()
        p = subprocess.run([str(executable),str(table),f.name],capture_output=True,
                           text=True,check=True)
    got = [float(l.split()[1]) for l in p.stdout.splitlines() if l.startswith('RESPONSE')]
    if len(got) != len(refs):
        raise ValueError('consumer did not return every sample')
    samples = list(zip(got,refs))
    scale = max(1e-30,max(abs(ref) for _,ref in samples))
    maxerr = max(abs(got-ref)/scale for got,ref in samples)
    out = dict(native_consumer_max_error=maxerr, tolerance=0.02, samples=len(samples),
               pitch_range=[float(min(pitches)),float(max(pitches))],
               bands='passing, trapped-passing boundary, incomplete turning band',
               normalization='max absolute native sample; finite-volume pitch-band means')
    if maxerr > out['tolerance']:
        raise ValueError(f'APPLETS versus native NEO-2 distribution: {out}')
    return out


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('exports',type=Path,nargs='+')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--table',type=Path)
    p.add_argument('--er-cgs',type=float)
    p.add_argument('--neo2-output',type=Path)
    p.add_argument('--n-theta',type=int,default=181)
    p.add_argument('--n-lambda',type=int,default=401)
    p.add_argument('--consumer',type=Path)
    a = p.parse_args()
    result = {str(f): oracles(read_fax(f)) for f in a.exports}
    er, epar = a.er_cgs, 0.
    if a.neo2_output:
        with h5py.File(a.neo2_output) as f:
            er = float(np.asarray(f['Er'][()]).item())
            if 'avEparB_ov_avb2' in f:
                epar = float(np.asarray(f['avEparB_ov_avb2'][()]).item())
    if a.neo2_output:
        result['dimensional_normalization'] = [dimensional_oracle(
            read_fax(f), a.neo2_output, result[str(f)]['qflux_reconstructed'])
            for f in a.exports]
    if a.table:
        if er is None:
            p.error('--table requires --er-cgs or --neo2-output; no implicit electric-field normalization')
        receipt = write_table(a.exports,a.table,er,epar,a.n_theta,a.n_lambda)
        result['table'] = receipt
        if a.consumer:
            result['consumer'] = [consumer_oracle(read_fax(f),a.table.resolve(),a.consumer.resolve(),
                                                   np.asarray(receipt['source_amplitudes']))
                                  for f in a.exports]
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({str(k): {n:v for n,v in r.items() if not n.startswith('qflux')}
                      if isinstance(r,dict) else r for k,r in result.items()},indent=2))


if __name__ == '__main__':
    main()
