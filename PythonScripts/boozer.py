#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 10:17:41 2016

@author: Christopher Albert
"""

def write_boozer_head(filename, version, shot, m0b, n0b, nsurf_fmt, nfp, psi_tor_a, aminor, Rmajor):
  import getpass

  global_variables_head_line=" m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]"

  with open(filename, 'w') as outfile:
    outfile.write("CC Boozer-coordinate data file\n")
    outfile.write("CC Version: " + version + '\n')
    outfile.write("CC Author: " + getpass.getuser() + '\n')
    outfile.write("CC shot: {:4d}".format(shot) + '\n')
    outfile.write(global_variables_head_line + '\n')
    outfile.write("{:4d}  {:4d}  {:3d} {:4d}  ".format(m0b, n0b, nsurf_fmt, nfp))
    outfile.write("{:12.6e}   {:12.6e}   {:12.6e}".format(psi_tor_a, aminor, Rmajor))
    outfile.write("\n")

def append_boozer_block_head(filename, s, iota, bsubvB, bsubuB, pprime, vp, enfp):
  import numpy as np
  import scipy.constants

  mu0 = scipy.constants.mu_0

  with open(filename, 'a') as outfile:
    outfile.write('        s               iota           Jpol/nper          '+
    'Itor            pprime         sqrt g(0,0)\n')
    outfile.write('                                          [A]           '+
    '[A]             [Pa]         (dV/ds)/nper\n')
    outfile.write(' {:16.8e}'.format(s))
    outfile.write(' {:16.8e}'.format(iota))
    outfile.write(' {:16.8e}'.format(-2.0*np.pi/mu0*bsubvB/enfp))
    outfile.write(' {:16.8e}'.format(-2.0*np.pi/mu0*bsubuB))
    outfile.write(' {:16.8e}'.format(pprime))
    outfile.write(' {:16.8e}'.format(-4.0*np.pi**2*vp/enfp))
    outfile.write('\n')

def append_boozer_block(filename, mb, nb, rmnb, zmnb, vmnb, bmnb, enfp):
  with open(filename, 'a') as f:
    f.write('    m    n      rmnc [m]         rmns [m]         zmnc [m]  '+
            '       zmns [m]         vmnc [m]         vmns [m]         '+
            'bmnc [T]         bmns [T]\n')
    for k in range(len(mb)):
      f.write(' {:4d} {:4d}'.format(mb[k],int(nb[k]/enfp)))
      f.write(' {:16.8e} {:16.8e}'.format(
          float((rmnb[k].real)),-float(rmnb[k].imag)))
      f.write(' {:16.8e} {:16.8e}'.format(
          float((zmnb[k].real)),-float(zmnb[k].imag)))
      f.write(' {:16.8e} {:16.8e}'.format(
          float((vmnb[k].real)),-float(vmnb[k].imag)))
      f.write(' {:16.8e} {:16.8e}'.format(
          float((bmnb[k].real)),-float(bmnb[k].imag)))
      f.write('\n')

def convert_to_boozer(infile, ks, outfile):
  # Import modules.

  # Look for a fourier transformation module.
  try:
    from fourier import fourierseries as fourierseries1
    print('Using f2py fourierseries .so')
  except:
    try:
      from fourier_win import fourierseries as fourierseries1
      print('Using Win32 f2py fourierseries DLL')
    except:
      def fourierseries1(fmn, u, v, m, n):
        ycpl = np.sum(fmn*np.exp(1j*(m*u + n*v)))
        return np.real(ycpl)
      print('Using Python fourierseries (SLOW!)')
  # Standard stuff.
  import math
  import numpy as np
  import scipy.interpolate as ip
  import string
  import sys
  import time
  #~ import scipy.constants.mu_0 as mu0

  def nextpow2(i):
    n = 1
    while n < i: n *= 2
    return n

  def fourierseries(fmn, u, v, m, n):
    if type(u) == float or u.size == 1:
      ycpl = np.sum(fmn*np.exp(1j*(m*u + n*v)))
      return np.real(ycpl)
    else:
      y = np.zeros(u.shape)
      for k in range(u.size):
        y[k] = fourierseries1(fmn, u[k], v[k], m, n)
      return y

  pi = np.pi

  #~ infile = sys.argv[1]
  #~ ks     = int(sys.argv[2])  # flux surface nr.
  nl     = 3          # number of points on each side for Lagrange interpolation
  plot   = False

  t = time.time()

  #infile = 'wout.test'

  low = True # use lower mode number for NetCDF (not Nyquist one)
  if infile.endswith('nc'):
    print('Reading NetCDF file {} at flux surface {}'.format(infile,ks))
    import scipy.io.netcdf as ncdf
    f = ncdf.netcdf_file(infile)
    data = f.variables

    enrho = np.copy(data['ns'].data)
    mlow = np.array(np.copy(data['xm'].data),int)
    nlow = np.array(-np.copy(data['xn'].data),int)
    m = np.array(np.copy(data['xm_nyq'].data),int)
    n = np.array(-np.copy(data['xn_nyq'].data),int)
    empol = int(np.max(np.abs(m)) + 1)
    entor = int(np.max(np.abs(n)))
    empoll = int(np.max(np.abs(mlow)) + 1)
    entorl = int(np.max(np.abs(nlow)))
    phip = -np.copy(data['phipf'].data)/(2.0*np.pi)
    pres = np.copy(data['presf'].data)
    empmnt= np.copy(data['mnmax_nyq'].data)
    empmntl= np.copy(data['mnmax'].data)
    buco =  np.copy(data['buco'].data)
    bvco =  np.copy(data['bvco'].data)
    iota =  np.copy(data['iotas'].data)
    enfp =  np.copy(data['nfp'].data)
    vp =    np.copy(data['vp'].data)

    rmnl = np.copy(data['rmnc'].data - 1.0j*data['rmns'].data)
    zmnl = np.copy(data['zmnc'].data - 1.0j*data['zmns'].data)
    lmnl = np.copy(data['lmnc'].data - 1.0j*data['lmns'].data)

    bsubumn = np.copy(data['bsubumnc'].data - 1.0j*data['bsubumns'].data)
    bsubvmn = np.copy(data['bsubvmnc'].data - 1.0j*data['bsubvmns'].data)
    bsubsmn = np.copy(data['bsubsmnc'].data - 1.0j*data['bsubsmns'].data)

    bsupumn = np.copy(data['bsupumnc'].data - 1.0j*data['bsupumns'].data)
    bsupvmn = np.copy(data['bsupvmnc'].data - 1.0j*data['bsupvmns'].data)

    del data
    f.close()
    # use only modes where all quantities are defined
    condi = (np.abs(m)<empoll) & (np.abs(n)<=entorl)

    if low:
      m = mlow
      n = nlow
      empol = empoll
      entor = entorl
      empmnt = empmntl
      rmn = rmnl
      zmn = zmnl
      lmn = lmnl
      bsubumn = bsubumn[:,condi]
      bsubvmn = bsubvmn[:,condi]
      bsubsmn = bsubsmn[:,condi]
      bsupumn = bsupumn[:,condi]
      bsupvmn = bsupvmn[:,condi]
    else:
      rmn = np.zeros(bsubumn.shape, complex)
      zmn = np.zeros(bsubumn.shape, complex)
      lmn = np.zeros(bsubumn.shape, complex)
      rmn[:,condi] = rmnl
      zmn[:,condi] = zmnl
      lmn[:,condi] = lmnl

  else:
    print('Reading text file {} at flux surface {}'.format(infile,ks))
    with open(infile) as f:
      lines = f.readlines()

    block = np.fromstring(lines[0], sep=' ')
    gamma = block[0]
    enfp  = int(block[1])
    enrho = int(block[2])
    block = np.fromstring(lines[1], sep=' ')
    empol = int(block[0])
    entor = int(block[1])
    empmnt= int(block[2])
    block = np.fromstring(lines[2], sep=' ')
    eiasym = int(block[0])
    phiedge  = block[1]

    lfourier = lines[3:enrho*empmnt+3]
    lprofile = lines[enrho*empmnt+3:enrho*empmnt+int(math.ceil((enrho-1)*12/5))+4]
    fourier = np.fromstring(string.join(lfourier),sep=' ').reshape(enrho,empmnt,16)
    profile = np.fromstring(string.join(lprofile), sep=' ')[0:12*(enrho-1)]\
        .reshape(enrho-1,12)

    # Fourier quantities
    rmn = fourier[:,:,0] - 1.0j*fourier[:,:,2]
    zmn = fourier[:,:,3] - 1.0j*fourier[:,:,1]
    bsupumn = fourier[:,:,4] - 1.0j*fourier[:,:,6]
    bsupvmn = fourier[:,:,5] - 1.0j*fourier[:,:,7]
    lmn = fourier[:,:,9] - 1j*fourier[:,:,8]
    bsubsmn = fourier[:,:,15] - 1.0j*fourier[:,:,12]
    bsubumn = fourier[:,:,10] - 1.0j*fourier[:,:,13]
    bsubvmn = fourier[:,:,11] - 1.0j*fourier[:,:,14]

    # Profile quantities
    iota = profile[:,0] # TODO: this seems to be full mesh, need Lagrange
    mass = profile[:,1]
    pres = profile[:,2]
    phip = profile[:,3]
    buco = profile[:,4]
    bvco = profile[:,5]
    phi  = profile[:,6]
    vp   = profile[:,7]
    overr = profile[:,8]
    jcuru = profile[:,9]
    jcurv = profile[:,10]
    specw = profile[:,11]

    # Modenumbers
    k = 0
    m = np.zeros(empmnt,dtype=int)
    n = np.zeros(empmnt,dtype=int)
    for mk in range(empol):
      nmin0 = -entor
      if mk < 1:
        nmin0 = 0
      for nk in range(nmin0,entor+1):
        m[k] = mk
        n[k] = -nk*enfp
        k = k+1

    # Adjustments
    iota = np.insert(iota, 0, 0.0)
    phip = np.insert(phip, 0, 0.0)
    buco = np.insert(buco, 0, 0.0)
    bvco = np.insert(bvco, 0, 0.0)
    vp = np.insert(vp, 0, 0.0)

  ns = enrho - 1
  ds = 1.0/ns
  s   = (np.arange(0,ns)+0.5)*ds
  sf  = (np.arange(0,ns+1))*ds

  if(ks < 3 or ks > ns-3):
    nl = 2

  if(ks < 2 or ks > ns-2):
    nl = 1

  print('Defining parameters and functions')
  #~ mu0 = 4e-7*np.pi

  cond1   = (m != 0)
  cond2   = (n != 0)
  m1 = m[cond1]; n1 = n[cond1]
  m2 = m[cond2]; n2 = n[cond2]

  s           = np.insert(s, 0, 0.0)
  #s = np.insert(s, len(s), 0.0)
  dpsitords   = -phip

  #%% Full mesh quantities
  ppoly  = ip.lagrange(sf[ks-nl:ks+nl],pres[ks-nl:ks+nl])
  pspoly = np.polyder(ppoly)
  psval  = np.polyval(ppoly,s[ks])

  rpoly = []; rmnval = []; rspoly = []; rsmnval = []
  zpoly = []; zmnval = []; zspoly = []; zsmnval = []
  lpoly = []; lmnval = []; lspoly = []; lsmnval = []
  for km in range(empmnt):
    rpoly.append(ip.lagrange(sf[ks-nl:ks+nl],rmn[ks-nl:ks+nl,km]))
    rspoly.append(np.polyder(rpoly[km]))
    rmnval.append(np.polyval(rpoly[km],s[ks]))
    rsmnval.append(np.polyval(rspoly[km],s[ks]))

    zpoly.append(ip.lagrange(sf[ks-nl:ks+nl],zmn[ks-nl:ks+nl,km]))
    zspoly.append(np.polyder(zpoly[km]))
    zmnval.append(np.polyval(zpoly[km],s[ks]))
    zsmnval.append(np.polyval(zspoly[km],s[ks]))

    lpoly.append(ip.lagrange(sf[ks-nl:ks+nl],lmn[ks-nl:ks+nl,km]))
    lspoly.append(np.polyder(lpoly[km]))
    lmnval.append(np.polyval(lpoly[km],s[ks]))
    lsmnval.append(np.polyval(lspoly[km],s[ks]))

  rmnval = np.array(rmnval); zmnval = np.array(zmnval); lmnval = np.array(lmnval)
  rsmnval= np.array(rsmnval);zsmnsval=np.array(zsmnval);lsmnval=np.array(lsmnval)

  # Cylindrical coordinates
  def r(u,v): return fourierseries(rmnval,u,v,m,n)
  def drdu(u,v): return fourierseries(1j*m*rmnval,u,v,m,n)
  def drds(u,v): return fourierseries(rsmnval,u,v,m,n)

  def z(u,v): return fourierseries(zmnval,u,v,m,n)
  def dzdu(u,v): return fourierseries(1j*m*zmnval,u,v,m,n)
  def dzds(u,v): return fourierseries(zsmnsval,u,v,m,n)

  # Stream function
  def lam(u,v): return fourierseries(lmnval,u,v,m,n)
  def dlamdu(u,v): return fourierseries(1j*m*lmnval,u,v,m,n)
  def dlamdv(u,v): return fourierseries(1j*n*lmnval,u,v,m,n)
  def dlamds(u,v): return fourierseries(lsmnval,u,v,m,n)

  def bsupu(u,v,fsi): return fourierseries(bsupumn[fsi,:],u,v,m,n)
  def bsupv(u,v,fsi): return fourierseries(bsupvmn[fsi,:],u,v,m,n)

  def bsubs(u,v,fsi): return fourierseries(bsubsmn[fsi,:],u,v,m,n)
  def bsubu(u,v,fsi): return fourierseries(bsubumn[fsi,:],u,v,m,n)
  def bsubv(u,v,fsi): return fourierseries(bsubvmn[fsi,:],u,v,m,n)

  def bmod2(u,v,fsi): return bsupu(u,v,fsi)*bsubu(u,v,fsi) +\
      bsupv(u,v,fsi)*bsubv(u,v,fsi)
  def bmod(u,v): return np.sqrt(bmod2(u,v,ks))

  # Metric tensor
  def G(u,v): return drdu(u,v)*dzds(u,v)-drds(u,v)*dzdu(u,v)
  def sqrtg(u,v): return np.abs(r(u,v)*G(u,v))

  # Alternative definition of stream function
  def dlamdu0(u,v): return sqrtg(u,v)/dpsitords[ks]*bsupv(u,v,ks)-1.0
  def dlamdv0(u,v): return (-sqrtg(u,v)/(iota[ks]*dpsitords[ks])*bsupu(u,v,ks)+1.0)*iota[ks]

  # redifine lambda
  #def dlamdu(u,v): return dlamdu0(u,v)
  #def dlamdv(u,v): return dlamdv0(u,v)

  # VMEC magnetic coordinates
  def uf(u,v,fsi): return u + lam(u,v)
  def sqrtgf(u,v): return sqrtg(u,v)/np.abs(1+dlamdu(u,v))
  def bsupuf(u,v,fsi): return (1+dlamdu(u,v,fsi))*bsupu(u,v,fsi)+\
      dlamdv(u,v,fsi)*bsubv(u,v,fsi)

  # %% Boozer coordinates
  bsubuB   = np.real(bsubumn[:,0])
  bsubvB   = np.real(bsubvmn[:,0])
  bsubuB2  = buco
  bsubvB2  = bvco
  pprime   = 0.0

  # output flux surface quantities
  # TODO: check sign convention. This one matches Strumberger output
  # TODO: check dV/ds and bvco enfp factor. Documentation says it's already included
  #       but it seems it needs to be added to get the correct output
  # TODO: pprime
  #~ with open('boozer{}.out'.format(ks), 'w') as outfile:
    #~ outfile.write('        s               iota           Jpol/nper          '+
    #~ 'Itor            pprime         sqrt g(0,0)\n')
    #~ outfile.write('                                          [A]           '+
    #~ '[A]             [Pa]         (dV/ds)/nper\n')
    #~ outfile.write(' {:16.8e}'.format(s[ks]))
    #~ outfile.write(' {:16.8e}'.format(iota[ks]))
    #~ outfile.write(' {:16.8e}'.format(-2.0*np.pi/mu0*bsubvB[ks]/enfp))
    #~ outfile.write(' {:16.8e}'.format(-2.0*np.pi/mu0*bsubuB[ks]))
    #~ outfile.write(' {:16.8e}'.format(pprime))
    #~ outfile.write(' {:16.8e}'.format(-4.0*np.pi**2*vp[ks]/enfp))
    #~ outfile.write('\n')
  append_boozer_block_head(outfile, s[ks], iota[ks], bsubvB[ks], bsubuB[ks], pprime, vp[ks], enfp)

  def sqrtgB(u,v): return np.abs(dpsitords[ks]*(iota[ks]*bsubuB[ks]+
                                  bsubvB[ks])/bmod2(u,v,ks))


  def hcheck(u,v): return sqrtgf(u,v)/sqrtgB(u,v) - 1.0

  # Boozer conversion
  print('Boozer conversion after Nuehrenberg/Zille')

  hmn1 = (bsubumn[ks,cond1]-1j*m1*bsubuB[ks]*lmnval[cond1])/\
          (1j*m1*(bsubuB[ks]+bsubvB[ks]/iota[ks]))
  hmn2 = (bsubvmn[ks,cond2]-1j*n2*bsubuB[ks]*lmnval[cond2])/\
          (1j*n2*(bsubuB[ks]+bsubvB[ks]/iota[ks]))

  hmn = np.zeros(m.shape, dtype='complex')
  hmn[cond2] = hmn2
  hmn[cond1] = hmn1

  def H(u,v): return fourierseries(hmn,u,v,m,n)
  def dHdu(u,v): return fourierseries(1j*m*hmn,u,v,m,n)
  def dHdv(u,v): return fourierseries(1j*n*hmn,u,v,m,n)

  def dth(u,v): return fourierseries(lmnval+hmn, u, v, m, n)
  def dph(u,v): return fourierseries(hmn/iota[ks], u, v, m, n)

  # Calculate Boozer modes
  m0b = 2*empol-1
  n0b = 4*entor+1
  mb = np.zeros((m0b-1)*n0b+(n0b-1)/2+1,dtype=int)
  nb = np.zeros((m0b-1)*n0b+(n0b-1)/2+1,dtype=int)
  k=0
  for mk in range(m0b):
    nmin0 = -2*entor
    if mk < 1:
      nmin0 = 0
    for nk in range(nmin0,2*entor+1):
      mb[k] = mk
      nb[k] = nk*enfp
      k = k+1

  nu = 6*np.max(np.abs(m))-1
  nv = 6*np.max(np.abs(n))-1
  du = 2.0*pi/nu
  dv = 2.0*pi/nv
  up = np.arange(0,2*pi,du)
  vp = np.arange(0,2*pi,dv)
  [U,V] = np.meshgrid(up,vp)
  U = U.flatten().T; V = V.flatten().T
  THB = U + dth(U,V)
  PHB = V + dph(U,V)
  R = r(U,V)
  Z = z(U,V)
  B = bmod(U,V)
  H1 = H(U,V)
  HCHECK = hcheck(U,V)
  Jbinv = np.abs((1.0+dlamdu(U,V)+dHdu(U,V))*(1.0+dHdv(U,V)/iota[ks])\
          -dHdu(U,V)/iota[ks]*(dlamdv(U,V)+dHdv(U,V)))

  print("Computing Boozer modes")
  dthmnb = np.zeros(mb.shape, complex)
  dphmnb = np.zeros(mb.shape, complex)
  rmnb = np.zeros(mb.shape, complex)
  zmnb = np.zeros(mb.shape, complex)
  bmnb = np.zeros(mb.shape, complex)
  hmnb = np.zeros(mb.shape, complex)
  hcheckmnb = np.zeros(mb.shape, complex)
  for km in range(len(mb)):
    efun = 2.0/(2.0*np.pi)**2*np.exp(-1.0j*(mb[km]*THB + nb[km]*PHB))*du*dv
    if (mb[km]==0) and (nb[km]==0):
      efun = efun/2.0
    dthmnb[km] = np.sum(efun*Jbinv*(THB-U))
    dphmnb[km] = np.sum(efun*Jbinv*(PHB-V))
    rmnb[km] = np.sum(efun*Jbinv*R)
    zmnb[km] = np.sum(efun*Jbinv*Z)
    bmnb[km] = np.sum(efun*Jbinv*B)
    hmnb[km] = np.sum(efun*Jbinv*H1)
    hcheckmnb[km] = np.sum(efun*Jbinv*HCHECK)

  vmnb = -enfp*dphmnb/(2*np.pi)

  append_boozer_block(outfile, mb, nb, rmnb, zmnb, vmnb, bmnb, enfp)
  #~ with open('boozer{}.out'.format(ks), 'a') as f:
    #~ f.write('    m    n      rmnc [m]         rmns [m]         zmnc [m]  '+
            #~ '       zmns [m]         vmnc [m]         vmns [m]         '+
            #~ 'bmnc [T]         bmns [T]\n')
    #~ for k in range(mb.size):
      #~ f.write(' {:4d} {:4d}'.format(mb[k],int(nb[k]/enfp)))
      #~ f.write(' {:16.8e} {:16.8e}'.format(
          #~ float((rmnb[k].real)),-float(rmnb[k].imag)))
      #~ f.write(' {:16.8e} {:16.8e}'.format(
          #~ float((zmnb[k].real)),-float(zmnb[k].imag)))
      #~ f.write(' {:16.8e} {:16.8e}'.format(
          #~ float((vmnb[k].real)),-float(vmnb[k].imag)))
      #~ f.write(' {:16.8e} {:16.8e}'.format(
          #~ float((bmnb[k].real)),-float(bmnb[k].imag)))
      #~ f.write('\n')

  elapsed = time.time() - t
  print('Elapsed time: {} s'.format(elapsed))


  #~ plot = True
  if plot:
    import matplotlib.pyplot as plt
    # Check contravariant B_Boozer
    u = np.linspace(-np.pi,np.pi,100)
    v = 0

    plt.figure(1)
    plt.clf()
    plt.subplot(1,2,1)
    plt.title('dlamdu')
    plt.plot(u, [dlamdu(ui,v) for ui in u])
    plt.plot(u, [dlamdu0(ui,v) for ui in u],'--')
    plt.legend(['dlamdu','dlamdu0'])

    plt.subplot(1,2,2)
    plt.title('dlamdv')
    plt.plot(u, [dlamdv(ui,v) for ui in u])
    plt.plot(u, [dlamdv0(ui,v) for ui in u],'--')
    plt.legend(['dlamdv','dlamdv0'])

    plt.figure(3)
    plt.clf()
    plt.title('sqrtg')
    plt.plot(u, [sqrtg(ui,v) for ui in u])
    plt.plot(u, [sqrtgf(ui,v) for ui in u],'--')
    plt.plot(u, [sqrtgB(ui,v) for ui in u],'-.')

    plt.figure(4)
    plt.clf()
    plt.subplot(1,2,1)
    plt.title('bsubuB')
    plt.plot(s, bsubuB)
    plt.plot(s, bsubuB2,'--')
    plt.subplot(1,2,2)
    plt.title('bsubvB')
    plt.plot(s, bsubvB)
    plt.plot(s, bsubvB2,'--')

  #  plt.figure(4)
  #  plt.subplot(1,2,1)
  #  plt.title('bsupu')
  #  plt.plot(s, bsupumn)
  #  plt.subplot(1,2,2)
  #  plt.title('bsupv')
  #  plt.plot(s, bsupvmn)
  #  plt.show()


  #  plt.figure(5)
  #  plt.subplot(1,2,1)
  #  plt.title('rmn')
  #  plt.plot(rmnval[(n==0)&(m<5)])
  #  plt.subplot(1,2,2)
  #  plt.title('zmn')
  #  plt.plot(zmnval[(n==0)&(m<5)])
  #  plt.show()


  #  plt.figure(4)
  #  plt.subplot(1,2,1)
  #  plt.title('bsupu')
  #  plt.plot(s, bsupumn)
  #  plt.subplot(1,2,2)
  #  plt.title('bsupv')
  #  plt.plot(s, bsupvmn)
  #  plt.show()

class BoozerFile:
  """Storing the information in a boozer file.

  This class is designed to store the information of a boozer file.

  TODO:
  So far, this can only handle boozer files with a single toroidal mode
  number.
  """
  comments = []
  m0b = 0.0
  n0b = 0.0
  nsurf = 0
  nper = 0
  flux = 0.0 #[Tm^2]
  a = 0.0
  R = 0.0 # [m]

  s = []
  iota = []
  Jpol_divided_by_nper = []
  Itor = []
  pprime = []
  sqrt_g_00 = []

  rmnc = [] # [m]
  rmns = [] # [m]
  zmnc = [] # [m]
  zmns = [] # [m]
  vmnc = []
  vmns = []
  bmnc = [] # [T]
  bmns = [] # [T]

  def read_boozer(self, filename: str):
    """Reads information from a file, whose name is given as a string.

    This routine will read the content of file 'filename', assuming it
    is a boozer file (thus expecting a specific format).
    The comments are available in a list, each line an entry.
    Global fields are available as simple elements.
    Fields that depend on radius are lists, those that depend on radius
    and mode number lists of lists.

    TODO:
    Be able to read in data with more than one 'n' mode number.
    """
    with open(filename) as f:
      lines = f.readlines()

    # Get header information, e.g. comments and sizes.
    for lineindex in range(len(lines)):
      line = lines[lineindex]
      if (len(line) > 2):
        # Comments start with 'CC' followed by a whitespace.
        if (line.split()[0] == 'CC'):
          self.comments.append(line[2+1:]) # 2+1 to ignore 'CC' and the following whitespace.

        if (lineindex > 1):
          if ((lines[lineindex-1].split())[0] == 'm0b'):
            self.m0b = int((lines[lineindex].split())[0])
            self.n0b = int((lines[lineindex].split())[1])
            self.nsurf = int((lines[lineindex].split())[2])
            self.nper = int((lines[lineindex].split())[3])
            self.flux = float((lines[lineindex].split())[4])
            self.a = float((lines[lineindex].split())[5])
            self.R = float((lines[lineindex].split())[6])
            break

    blocklines = []

    blockindex = -1
    for lineindex in range(len(lines)):
      line = lines[lineindex]
      if (line.split()[0] == 's'):
        blockindex += 1
        blocklines.append([])

      if (blockindex >= 0 and len(line) > 0):
        blocklines[blockindex].append(line)

    if (len(blocklines) != self.nsurf):
      print('m0b = ' + str(self.m0b))
      print('n0b = ' + str(self.n0b))
      print('nsurf = ' + str(self.nsurf))
      print('nper = ' + str(self.nper))
      print('flux = ' + str(self.flux))
      print('a = ' + str(self.a))
      print('R = ' + str(self.R))
      print(str(len(blocklines)) + ' != ' + str(self.nsurf))
      raise Exception

    head_number_of_lines = 4

    for i in range(self.nsurf):
      if (len(blocklines[i]) != self.m0b + 1 + head_number_of_lines): # +1 for zero mode
        raise Exception

      self.rmnc.append([])
      self.rmns.append([])
      self.zmnc.append([])
      self.zmns.append([])
      self.vmnc.append([])
      self.vmns.append([])
      self.bmnc.append([])
      self.bmns.append([])

      for j in range(head_number_of_lines, self.m0b + 1 + head_number_of_lines):
        if (j == 2):
          line_split = blocklines[i][j].split();
          self.s.append(float(line_split[0]))
          self.iota.append(float(line_split[1]))
          self.Jpol_divided_by_nper.append(float(line_split[2]))
          self.Itor.append(float(line_split[3]))
          self.pprime.append(float(line_split[4]))
          self.sqrt_g_00.append(float(line_split[5]))

        if (j > 3):
          line_split = blocklines[i][j].split();
          self.rmnc[i].append(float(line_split[2]))
          self.rmns[i].append(float(line_split[3]))
          self.zmnc[i].append(float(line_split[4]))
          self.zmns[i].append(float(line_split[5]))
          self.vmnc[i].append(float(line_split[6]))
          self.vmns[i].append(float(line_split[7]))
          self.bmnc[i].append(float(line_split[8]))
          self.bmns[i].append(float(line_split[9]))

  def __init__(self, filename: str):
    """Init routine which takes a string, representing the file to read.
    """
    self.read_boozer(filename)

if __name__ == "__main__":
  import sys

  import getHeadDataVMEC

  if (len(sys.argv) < 2):
    print("Usage:")
    print("./boozer.py infilename fluxsurfacenumber")
  else:
    infile = sys.argv[1]
    nsurf = int(sys.argv[2])
    wout_name = infile + '.bc'

    [nfp, psi_tor_a, aminor, Rmajor, m0b, n0b] = getHeadDataVmecNc(infile)

    write_boozer_head(wout_name, '01', m0b, n0b, nsurf, nfp, psi_tor_a, aminor, Rmajor)

    for ind in range(1, nsurf+1):
      convert_to_boozer(infile, ind, wout_name)
