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
%   gridpoints: number of points to use for the (equidistant) rho_pol grid.
%   do_plots: logical, if true, do some plots.
%   input_unit_type: determines in what units the input is.
%     1: 10^19/m^3 for density (=10^13/cm^3), keV for temperatures,
%       km?/s for velocity.
%     2: 1/m^3 for density, eV for temperatures, rad/s for velocity
%
% Output:
%  rho_pol
%  rho_tor
%  ne_si
%  Ti_eV
%  Te_eV
%  vrot
function [rho_pol, rho_tor, ne_si, Ti_eV, Te_eV, vrot] = load_profile_data(path_to_shot, data_source, gridpoints, do_plots, input_unit_type)


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
  otherwise
    transform_density = transform_1013overcm3_to_oneoverm3;
    transform_temperature = transform_keV_to_eV;
    transform_rotation = transform_kHz_to_Hz;
  end

  rho_pol = linspace(0,1,gridpoints);

  frp = load([path_to_shot, data_source.rhopoloidal.filename]);
  frt = load([path_to_shot, data_source.rhotoroidal.filename]);

  rho_tor = spline(frp(:, data_source.rhopoloidal.column), frt(:,data_source.rhotoroidal.column), rho_pol);
  rho_tor(1) = 0;
  rho_tor(end) = 1;

  frp = load([path_to_shot, data_source.electron_density.filename]);
  ne_si = spline(frp(:,1), frp(:, data_source.electron_density.column), rho_pol)*transform_density;


  frp = load([path_to_shot, data_source.electron_temperature.filename]);
  Te_eV = spline(frp(:,1), frp(:, data_source.electron_temperature.column), rho_pol)*transform_temperature;


  frp=load([path_to_shot, data_source.ion_temperature.filename]);
  rho_fit2=rho_pol.^2;

  % octave does not have a 'fit' function, only 'fsolve'. Thus this was replaced.
  %[fitobject,gof] = fit(frp(:,1).^2, frp(:,2) ,'poly6');
  %fit2 = fitobject.p1.*(rho_fit2.^6) + fitobject.p2.*(rho_fit2.^5) + fitobject.p3.*(rho_fit2.^4) + fitobject.p4.*(rho_fit2.^3) + fitobject.p5.*(rho_fit2.^2) + fitobject.p6.*rho_fit2+fitobject.p7;
  [fitobject, gof] = polyfit(frp(:,1).^2, frp(:, data_source.electron_temperature.column), 6);
  fit2 = polyval(fitobject, rho_fit2);

  Ti_eV = fit2*transform_temperature;

  if do_plots
    figure
    plot(rho_pol,Ti_eV*transform_eV_to_keV,'r',frp(:,1),frp(:,3),'b--',frp(:,1),frp(:,4),'b--')
    xlabel('\rho_{pol}')
    ylabel('T_i  [keV]')
  end


  frp = load([path_to_shot, data_source.rotation_velocity.filename]);
  %frr = load([path_to_shot, data_source.major_radius.filename]);

  switch (input_unit_type)
  case 1
    frp(:, data_source.rotation_velocity.column) = frp(:, data_source.rotation_velocity.column)./frr(:, data_source.major_radius.column);
  case 2
    frp(:, data_source.rotation_velocity.column) = frp(:, data_source.rotation_velocity.column);
  otherwise
    frp(:, data_source.rotation_velocity.column) = frp(:, data_source.rotation_velocity.column)./frr(:, data_source.major_radius.column);
  end

  frp(:, data_source.rotation_velocity.column) = frp(:, data_source.rotation_velocity.column)*transform_rotation;

  [fitobject, gof] = polyfit(frp(:,1).^2, frp(:, data_source.rotation_velocity.column), 6);
  fit2 = polyval(fitobject, rho_fit2);

  vrot = fit2;

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

end
