def get_D_coef_neo2par(neo2par_runfolder, T, m, Z):
    import os
    import numpy as np
    import h5py

    c = 3e10 # cm/s
    e = -4.8e-10 # statC
    transp_par = os.path.join(neo2par_runfolder, "fulltransp.h5")
    magnetics_par = os.path.join(neo2par_runfolder, "efinal.h5")

    with h5py.File(transp_par) as file:
        qflux = file["qflux"][()]
        int_1_ov_B_dl = file["dl1obhat"][()]
        average_nabla_s = file["av_nabla_stor"][()]
    with h5py.File(magnetics_par) as file:
        Bref = file["bmod0"][()] * 1e4 # G

    v_T = np.sqrt(2 * T / m)
    omega_c = np.abs(Z) * e * Bref / (m * c)
    rho = v_T/omega_c
    beta = np.array([rho/average_nabla_s, rho/average_nabla_s, Bref])

    Dij = np.zeros_like(qflux)

    remap = [0, 2, 1] # qflux has switched heat flux/bootstrap current index
    for i in range(3):
        for j in range(3):
            Dij[i,j] = - beta[i] * beta[j] * qflux[remap[i],remap[j]] * v_T / int_1_ov_B_dl

    return Dij