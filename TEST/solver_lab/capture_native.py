#!/usr/bin/env python3
"""Build-only instrumentation: capture native QL algebra, not a new discretization.

Only public example geometry and documented namelists are used. NEO2_CAPTURE=AX
or NA exits collectively before the selected sparse factorization. The working
Fortran source is modified for the experiment; use a clean checkout each time.
"""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path
import re
import shutil

HOOK = r'''
    subroutine capture_native_operator()
        use, intrinsic :: iso_fortran_env, only: int32
        character(len=32) :: wanted
        character(len=128) :: filename
        character(len=2) :: kind_name
        integer :: capture_unit, capture_status, j, nnz_original
        integer(int32) :: header(14)
        complex(dp), allocatable :: probe(:), ci_probe(:), original_probe(:)
        call get_environment_variable('NEO2_CAPTURE', wanted, status=capture_status)
        if (capture_status /= 0) return
        kind_name = 'NA'
        if (problem_type) kind_name = 'AX'
        if (trim(wanted) /= kind_name) return
        if (storage_size(1) /= 32 .or. storage_size(1.0_dp) /= 64) &
            error stop 'unsupported capture numeric kinds'
        nnz_original = nz
        if (problem_type) nnz_original = nz_symm
        header = int([1, n_2d_size, nz, nnz_original, lag+1, leg+1, &
            iend-ibeg+1, npart+1, num_spec, ispec, isw_intp, &
            merge(1,0,problem_type), merge(1,0,addboucol), 4], int32)
        write(filename, '(a,a,a,i0,a)') 'native_', kind_name, '_', ispec, '.bin'
        open(newunit=capture_unit, file=trim(filename), access='stream', &
             form='unformatted', status='replace', action='write')
        write(capture_unit) 'NEO2CAP1', header
        write(capture_unit) int(irow(1:nz),int32), int(icol(1:nz),int32), amat_sp(1:nz)
        if (problem_type) then
            write(capture_unit) int(irow_symm,int32), int(icol_symm,int32), &
                cmplx(amat_symm,0.0_dp,dp)
        else
            write(capture_unit) int(irow(1:nz),int32), int(icol(1:nz),int32), amat_sp(1:nz)
        end if
        write(capture_unit) source_vector_all, flux_vector
        write(capture_unit) int(npl,int32), int(ind_start,int32), fact_pos_e, fact_neg_e
        write(capture_unit) pleg_bra(0:leg,:,:), pleg_ket(0:leg,:,:), ailmm_aa
        write(capture_unit) densvec_bra, densvec_ket, energvec_bra, energvec_ket
        allocate(probe(n_2d_size),ci_probe(n_2d_size),original_probe(n_2d_size))
        do j=1,n_2d_size
            probe(j)=cmplx(sin(real(j+17*ispec,dp)*0.317_dp), &
                          cos(real(3*j+ispec,dp)*0.193_dp),dp)
        end do
        ci_probe=(0.0_dp,0.0_dp)
        if (isw_intp == 1) call integral_part(probe,ci_probe)
        original_probe=-ci_probe
        if (problem_type) then
            do j=1,nz_symm
                original_probe(irow_symm(j))=original_probe(irow_symm(j)) &
                    + amat_symm(j)*probe(icol_symm(j))
            end do
        else
            do j=1,nz
                original_probe(irow(j))=original_probe(irow(j)) &
                    + amat_sp(j)*probe(icol(j))
            end do
        end if
        write(capture_unit) probe,ci_probe,original_probe
        close(capture_unit)
        print *, 'NEO2_NATIVE_CAPTURE',kind_name,ispec,header
        call mpro%barrier()
        call mpro%deinit()
        stop
    end subroutine capture_native_operator
'''

def instrument(root: Path) -> None:
    path=root/'NEO-2-QL/ripple_solver_ArnoldiOrder2_test.f90'
    text=path.read_text()
    if 'subroutine capture_native_operator()' in text:
        raise RuntimeError('already instrumented; use a clean checkout')
    anchor='        call column_full2pointer(icol(1:nz), ipcol)'
    anchor2='    subroutine integral_part(vec_in, vec_out)'
    if text.count(anchor)!=1 or text.count(anchor2)!=1:
        raise RuntimeError('native solve changed: review hook placement')
    text=text.replace(anchor,'        call capture_native_operator()\n\n'+anchor)
    text=text.replace(anchor2,HOOK+'\n'+anchor2)
    path.write_text(text)
    print('instrumented_sha256',hashlib.sha256(text.encode()).hexdigest())

def prepare(root: Path,dest: Path,eta: int,steps: int,lag: int,density: float) -> None:
    dest.mkdir(parents=True,exist_ok=True)
    for name in ('example_bg_field.bc','example_pert_field.bc'):
        shutil.copyfile(root/'DOC'/name,dest/name)
    shutil.copyfile(root/'DOC/neo.in.ql-full',dest/'neo.in')
    text=(root/'DOC/neo2.in.ql-full').read_text()
    replacements={
        'eta_part':str(eta),'eta_part_trapped':str(max(4,eta//2)),
        'nstep':str(steps),'nperiod':'1','lag':str(lag),
        'leg':'3','legmax':'5','omp_num_threads':'1',
        'in_file_pert':"'example_pert_field.bc'",'isw_calc_er':'0',
        'isw_calc_magdrift':'0','isw_mag_shear':'0','isw_ntv_mode':'0',
        'prop_write':'0','prop_reconstruct':'0','bsfunc_max_index':'3000',
        'isw_qflux_NA':'1','lsw_multispecies':'.true.','num_spec':'2',
        'conl_over_mfp':'1.0d-3','isw_coul_log':'1',
        'n_vec':f'2*{density:.16e}','t_vec':'2*1.602176634d-9',
        'collop_base_exp':'0','collop_base_prj':'0',
        'lsw_read_precom':'.false.','lsw_write_precom':'.false.',
    }
    for key,value in replacements.items():
        pattern=rf'(?im)^(\s*{re.escape(key)}\s*=)[^!\n]*'
        text,count=re.subn(pattern,lambda m:m[1]+' '+value+' ',text)
        if count!=1: raise RuntimeError(f'expected exactly one {key}, got {count}')
    (dest/'neo2.in').write_text(text)
    (dest/'capture-input-sha256.txt').write_text(hashlib.sha256(text.encode()).hexdigest()+'\n')

def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    ap.add_argument('--instrument',action='store_true')
    ap.add_argument('--prepare',type=Path)
    ap.add_argument('--eta',type=int,default=8)
    ap.add_argument('--steps',type=int,default=12)
    ap.add_argument('--lag',type=int,default=3)
    ap.add_argument('--density',type=float,default=1e13)
    a=ap.parse_args()
    if a.instrument: instrument(a.root)
    if a.prepare: prepare(a.root,a.prepare,a.eta,a.steps,a.lag,a.density)
if __name__=='__main__': main()
