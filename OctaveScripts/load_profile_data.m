% function [rho_pol, rho_tor, ne_si, Ti_eV, vrot] = load_profile_data(path_to_shot, shot_designation, do_plots)
%
% Function that loads profile data for a shot and makes some
% preprocessing.
%
% It requires the quantities rhotoroidal, electron density, electron
% temperature, ion temperature, rotation velocity and major radius, each
% as a function of rhopoloidal.
% Rhotoroidal might be given implicitly by having one quantity (e.g.
% electron density) as a function of both, rhotoroidal and rhopoloidal
% (maybe in different files), or explicitly in an own file.
% Major radius can either be an additional column in the file with the
% rotation velocity, or it can be given in an additional file with the
% same grid as the rotation velocity.
% In all the other cases it is assumed that the first column of the
% respective file is rhopoloidal.
%
% Input:
%   path_to_shot: string with the location to the shot, e.g. "example_folder/SHOTS/NUMBER_12345/".
%   data_source: a structure with at least the fields 'rhopoloidal',
%     'rhotoroidal', 'electron_density', 'electron_temperature',
%     'ion_temperature', 'rotation_velocity' and 'major_radius'
%     corresponding to the required quantities to generate the output.
%     Each of these fields is expected to be a structure, with at least
%     the fields 'filename' and 'column'. The former is a string that
%     contains the name off the file where to find the data for the
%     coresponding quantity, while the later is the number of the
%     column in which to find this quantity in the file.
%   gridpoints: either number of points to use for the (equidistant)
%     rho_pol grid (scalar), or array with three values, first number of
%     grid points, second lower boundary for flux surface label, and
%     third upper boundary.
%   do_plots: logical, if true, do some plots.
%   input_unit_type: determines in what units the input is.
%     1: 10^19/m^3 for density (=10^13/cm^3), keV for temperatures,
%       km?/s for velocity. (default)
%     2: 1/m^3 for density, eV for temperatures, rad/s for velocity
%   switch_grid: determines what grid to use.
%     1: equal spaced in rho poloidal (default)
%     2: equal spaced in rho toroidal
%     3: equal spaced in sqrt(rho poloidal), this means there will be
%       more points close to the edge, compared to the other two
%       options.
%
% Output:
%  rho_pol: rho_poloidal grid.
%  rho_tor: rho_toroidal at rho_pol.
%  ne_si: density at rho_pol.
%  Ti_eV: ion temperature at rho_pol.
%  Te_eV: electron temperature at rho_pol.
%  vrot: rotation velocitz at rho_pol.
function [rho_pol, rho_tor, ne_si, Ti_eV, Te_eV, vrot] = load_profile_data(path_to_shot, data_source, gridpoints, do_plots, input_unit_type, switch_grid)
  if nargin() < 5 || isempty(input_unit_type)
    input_unit_type = 1;
  end
  if nargin() < 6 || isempty(switch_grid)
    switch_grid = 1;
  end
  if (size(gridpoints(:), 1) == 3)
    number_gridpoints = gridpoints(1);
    lower_limit_flux_label = gridpoints(2);
    upper_limit_flux_label = gridpoints(3);
  else
    number_gridpoints = gridpoints;
    lower_limit_flux_label = 0.0;
    upper_limit_flux_label = 1.0;
  end

  write_data = 1;

  % Definitions for constants used for conversion
  transform_keV_to_eV = 1.0e3;
  transform_eV_to_keV = 1.0/transform_keV_to_eV;

  transform_kHz_to_Hz = 1.0e3;
  transform_Hz_to_kHz = 1.0/transform_kHz_to_Hz;

  transform_oneoverm3_to_1013overcm3 = 1.0e-19;
  transform_1013overcm3_to_oneoverm3 = 1.0/transform_oneoverm3_to_1013overcm3;

  switch (input_unit_type)
  case 1
    transform_density = transform_1013overcm3_to_oneoverm3;
    transform_temperature = transform_keV_to_eV;
    transform_rotation = transform_kHz_to_Hz;
  case 2
    transform_density = 1;
    transform_temperature = 1;
    transform_rotation = 1;
  end

  frp = load([path_to_shot, data_source.rhopoloidal.filename]);
  frt = load([path_to_shot, data_source.rhotoroidal.filename]);

  switch (switch_grid)
  case 1
    rho_pol = linspace(lower_limit_flux_label, upper_limit_flux_label, number_gridpoints);

    L = frp(:, data_source.rhopoloidal.column) <= 1.0;
    rho_tor = local_used_fit_function(frp(L, data_source.rhopoloidal.column), frt(L,data_source.rhotoroidal.column), rho_pol, 5);
    if (lower_limit_flux_label <= 0.0)
      rho_tor(1) = 0;
    end
    if (upper_limit_flux_label >= 1.0)
      rho_tor(end) = 1;
    end
  case 2
    rho_tor = linspace(0,1,number_gridpoints);

    rho_pol = local_used_fit_function(frt(:,data_source.rhotoroidal.column), frp(:, data_source.rhopoloidal.column), rho_tor, 5);
    rho_pol(1) = 0;
    rho_pol(end) = 1;
  case 3
    rho_pol = sqrt(linspace(lower_limit_flux_label, upper_limit_flux_label, number_gridpoints));

    L = frp(:, data_source.rhopoloidal.column) <= 1.0;
    rho_tor = local_used_fit_function(frp(L, data_source.rhopoloidal.column), frt(L,data_source.rhotoroidal.column), rho_pol, 5);
    if (lower_limit_flux_label <= 0.0)
      rho_tor(1) = 0;
    end
    if (upper_limit_flux_label >= 1.0)
      rho_tor(end) = 1;
    end
  end
  interpolation_grid = rho_pol.^2;

  frp = load([path_to_shot, data_source.electron_density.filename]);
  ne_si = local_used_fit_function(frp(:,1).^2, frp(:, data_source.electron_density.column), interpolation_grid, 6)*transform_density;


  frp = load([path_to_shot, data_source.electron_temperature.filename]);
  Te_eV = local_used_fit_function(frp(:,1).^2, frp(:, data_source.electron_temperature.column), interpolation_grid, 6)*transform_temperature;


  frp=load([path_to_shot, data_source.ion_temperature.filename]);
  Ti_eV = local_used_fit_function(frp(:,1).^2, frp(:, data_source.ion_temperature.column), interpolation_grid, 6)*transform_temperature;

  if do_plots
    figure
    plot(rho_pol,Ti_eV*transform_eV_to_keV,'r',frp(:,1),frp(:,3),'b--',frp(:,1),frp(:,4),'b--')
    xlabel('\rho_{pol}')
    ylabel('T_i  [keV]')
  end


  frp = load([path_to_shot, data_source.rotation_velocity.filename]);
  if (input_unit_type == 1)
    frr = load([path_to_shot, data_source.major_radius.filename]);
  end

  switch (input_unit_type)
  case 1
    frp(:, data_source.rotation_velocity.column) = frp(:, data_source.rotation_velocity.column)./frr(:, data_source.major_radius.column);
  case 2
    frp(:, data_source.rotation_velocity.column) = frp(:, data_source.rotation_velocity.column);
  otherwise
    frp(:, data_source.rotation_velocity.column) = frp(:, data_source.rotation_velocity.column)./frr(:, data_source.major_radius.column);
  end

  frp(:, data_source.rotation_velocity.column) = frp(:, data_source.rotation_velocity.column)*transform_rotation;

  vrot = local_used_fit_function(frp(:,1).^2, frp(:, data_source.rotation_velocity.column), interpolation_grid, 6);

  if do_plots
    figure
    plot(rho_pol,vrot/pi*transform_Hz_to_kHz,'r',frp(:,1),frp(:,3)/pi*1e-3,'b--',frp(:,1),frp(:,4)/pi*1e-3,'b--')
    title('f_{probe} = n\Omega_{tor}/(2 \pi),   n=2')
    xlabel('\rho_{pol}')
    ylabel('f_{probe}  [kHz]')

    figure
    plot(rho_pol,Te_eV*transform_eV_to_keV)
    xlabel('\rho_{pol}')
    ylabel('T_e  [keV]')

    figure
    plot(rho_pol,ne_si*transform_oneoverm3_to_1013overcm3)
    xlabel('\rho_{pol}')
    ylabel('n_e  [10^{13} cm^{-3}]')
  end

  if write_data
    fid = 0;
    fid = fopen([path_to_shot, 'flux_coordinates_densities_temperatures.dat'], 'w');
    fprintf(fid, "!              s          rho_tor          rho_pol density[1/cm^3] electrons    ions temperature[eV] electrons    ions\n");
    for k = 1:number_gridpoints
      fprintf(fid, "%16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e\n", rho_tor(k).^2, rho_tor(k), rho_pol(k), ne_si(k)*1.0e-6, ne_si(k)*1.0e-6, Te_eV(k), Ti_eV(k));
    end
    fclose(fid);
  end
end

% Local wrapper for fit functions
%
% A simple wrapper, to allow easy switching between different
% (polynomial, spline) interpolation methods.
%
% input:
% ------
% old_x: old grid on which the data is given.
% old_y: data values on old grid.
% new_x: new grid on which to calculate the y values.
% order: optional, order of polynomial interpolation (if used). [6]
% type: optional, type of interpolation to use. ['s']
%   's': spline
%   'p': polynomial
%   'b': both, blend smoothly from polynomial at the axis to spline at
%     the edge.
function fit = local_used_fit_function(old_x, old_y, new_x, order, type)
  if nargin() < 4 || isempty(order)
    order = 6;
  end
  if nargin() < 5 || isempty(type)
    type = 's';
  end

  switch (type)
  case 's'
    fit = spline(old_x, old_y, new_x);
  case 'p'
    fitobject = polyfit(old_x, old_y, order);
    fit = polyval(fitobject, new_x);
  case 'b'
    fitobject = polyfit(old_x, old_y, order);
    weight = (new_x - new_x(1)) ./ (new_x(end) - new_x(1));
    fit = weight .* spline(old_x, old_y, new_x) + (1 - weight) .* polyval(fitobject, new_x);
  otherwise
    print(['Error: unknown value for type of interpolation: ', type])
    quit()
  end
end
