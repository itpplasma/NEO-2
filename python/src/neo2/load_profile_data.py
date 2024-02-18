import numpy as np
import os

def local_used_fit_function(old_x, old_y, new_x, order=6, type='s'):
    if type == 's':
        fit = np.interp(new_x, old_x, old_y)
    elif type == 'p':
        fit = np.polyval(np.polyfit(old_x, old_y, order), new_x)
    elif type == 'b':
        fitobject = np.polyfit(old_x, old_y, order)
        weight = (new_x - new_x[0]) / (new_x[-1] - new_x[0])
        fit = weight * np.interp(new_x, old_x, old_y) + (1 - weight) * np.polyval(fitobject, new_x)
    else:
        print('Error: unknown value for type of interpolation:', type)
        quit()
    return fit

def load_profile_data(path_to_shot, data_source, gridpoints, do_plots=False, input_unit_type=1, switch_grid=1):
    if input_unit_type is None:
        input_unit_type = 1
    if switch_grid is None:
        switch_grid = 1

    if len(np.array(gridpoints).flatten()) == 3:
        number_gridpoints = gridpoints[0]
        lower_limit_flux_label = gridpoints[1]
        upper_limit_flux_label = gridpoints[2]
    else:
        number_gridpoints = gridpoints
        lower_limit_flux_label = 0.0
        upper_limit_flux_label = 1.0

    write_data = 1

    transform_keV_to_eV = 1.0e3
    transform_eV_to_keV = 1.0 / transform_keV_to_eV

    transform_krads_to_rads = 1.0e3
    transform_rads_to_krads = 1.0 / transform_krads_to_rads

    transform_oneoverm3_to_1013overcm3 = 1.0e-19
    transform_1013overcm3_to_oneoverm3 = 1.0 / transform_oneoverm3_to_1013overcm3

    if input_unit_type == 1:
        transform_density = transform_1013overcm3_to_oneoverm3
        transform_temperature = transform_keV_to_eV
        transform_rotation = transform_krads_to_rads
    elif input_unit_type == 2:
        transform_density = 1
        transform_temperature = 1
        transform_rotation = 1

    frp = np.loadtxt(os.path.join(path_to_shot, data_source['rhopoloidal']['filename']))
    frt = np.loadtxt(os.path.join(path_to_shot, data_source['rhotoroidal']['filename']))

    if switch_grid == 1:
        rho_pol = np.linspace(lower_limit_flux_label, upper_limit_flux_label, number_gridpoints)

        L = frp[:, data_source['rhopoloidal']['column']] <= 1.0
        rho_tor = local_used_fit_function(frp[L, data_source['rhopoloidal']['column']],
                                           frt[L, data_source['rhotoroidal']['column']], rho_pol, 5)
        if lower_limit_flux_label <= 0.0:
            rho_tor[0] = 0
        if upper_limit_flux_label >= 1.0:
            rho_tor[-1] = 1
    elif switch_grid == 2:
        rho_tor = np.linspace(0, 1, number_gridpoints)

        rho_pol = local_used_fit_function(frt[:, data_source['rhotoroidal']['column']],
                                          frp[:, data_source['rhopoloidal']['column']], rho_tor, 5)
        rho_pol[0] = 0
        rho_pol[-1] = 1
    elif switch_grid == 3:
        rho_pol = np.sqrt(np.linspace(lower_limit_flux_label, upper_limit_flux_label, number_gridpoints))

        L = frp[:, data_source['rhopoloidal']['column']] <= 1.0
        rho_tor = local_used_fit_function(frp[L, data_source['rhopoloidal']['column']],
                                           frt[L, data_source['rhotoroidal']['column']], rho_pol, 5)
        if lower_limit_flux_label <= 0.0:
            rho_tor[0] = 0
        if upper_limit_flux_label >= 1.0:
            rho_tor[-1] = 1

    interpolation_grid = rho_pol ** 2

    frp = np.loadtxt(os.path.join(path_to_shot, data_source['electron_density']['filename']))
    ne_si = local_used_fit_function(frp[:, 0] ** 2, frp[:, data_source['electron_density']['column']],
                                    interpolation_grid, 6) * transform_density

    frp = np.loadtxt(os.path.join(path_to_shot, data_source['electron_temperature']['filename']))
    Te_eV = local_used_fit_function(frp[:, 0] ** 2, frp[:, data_source['electron_temperature']['column']],
                                    interpolation_grid, 6) * transform_temperature

    frp = np.loadtxt(os.path.join(path_to_shot, data_source['ion_temperature']['filename']))
    Ti_eV = local_used_fit_function(frp[:, 0] ** 2, frp[:, data_source['ion_temperature']['column']],
                                    interpolation_grid, 6) * transform_temperature

    frp = np.loadtxt(os.path.join(path_to_shot, data_source['rotation_velocity']['filename']))
    if input_unit_type == 1:
        frr = np.loadtxt(os.path.join(path_to_shot, data_source['major_radius']['filename']))
    if input_unit_type == 1:
        frp[:, data_source['rotation_velocity']['column']] /= frr[:, data_source['major_radius']['column']]
    elif input_unit_type == 3:
        frp[:, data_source['rotation_velocity']['column']] /= frr[:, data_source['major_radius']['column']]
    frp[:, data_source['rotation_velocity']['column']] *= transform_rotation

    vrot = local_used_fit_function(frp[:, 0] ** 2, frp[:, data_source['rotation_velocity']['column']],
                                   interpolation_grid, 6)

    if write_data:
        with open(os.path.join(path_to_shot, 'flux_coordinates_densities_temperatures.dat'), 'w') as fid:
            fid.write("!              s          rho_tor          rho_pol density[1/cm^3] electrons    ions temperature[eV] electrons    ions\n")
            for k in range(number_gridpoints):
                fid.write(f"{rho_tor[k]**2:16.10e} {rho_tor[k]:16.10e} {rho_pol[k]:16.10e} "
                          f"{ne_si[k]*1.0e-6:16.10e} {ne_si[k]*1.0e-6:16.10e} {Te_eV[k]:16.10e} {Ti_eV[k]:16.10e}\n")

    return rho_pol, rho_tor, ne_si, Ti_eV, Te_eV, vrot