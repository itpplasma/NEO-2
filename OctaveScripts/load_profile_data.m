% function [rho_pol, rho_tor, ne_si, Ti_eV, vrot] = load_profile_data(path_to_shot, shot_designation, do_plots)
%
% Function that loads profile data for a shot and makes some
% preprocessing.
%
% It requires the six files ne_ida_shot_designation_rhopol.dat,
% ne_ida_shot_designation_rhotor.dat, Te_ida_shot_designation_rhopol.dat,
% Ti_cez_shot_designation_rhopol.dat, vrot_ida_shot_designation_rhopol.dat,
% and vrot_ida_shot_designation_R.dat. Here shot_designation is a common
% part in the filenames, that might be specific to the shot (and
% timepoint) one whant to load.
%
%
% Input:
%   path_to_shot: string with the location to the shot, e.g. "example_folder/SHOTS/NUMBER_12345/".
%   shot_designation: string whith the common part of the filenames. See
%     above for the expected naming scheme.
%   gridpoints: number of points to use for the (equidistant) rho_pol grid.
%   do_plots: logical, if true, do some plots.
%
% Output:
%  rho_pol
%  rho_tor
%  ne_si
%  Ti_eV
%  Te_eV
%  vrot
function [rho_pol, rho_tor, ne_si, Ti_eV, Te_eV, vrot] = load_profile_data(path_to_shot, shot_designation, gridpoints, do_plots)
  % Definitions for constants used for conversion
  transform_keV_to_eV = 1.0e3;
  transform_eV_to_keV = 1.0/transform_keV_to_eV;

  transform_kHz_to_Hz = 1.0e3;
  transform_Hz_to_kHz = 1.0/transform_kHz_to_Hz;

  transform_oneoverm3_to_1013overcm3 = 1.0e-19;
  transform_1013overcm3_to_oneoverm3 = 1.0/transform_oneoverm3_to_1013overcm3;

  rho_pol=linspace(0,1,gridpoints);

  frp=load([path_to_shot, '/ne_ida_', shot_designation,'_rhopol.dat']);
  frt=load([path_to_shot, '/ne_ida_', shot_designation,'_rhotor.dat']);

  rho_tor=spline(frp(:,1),frt(:,1),rho_pol);
  rho_tor(1)=0;
  rho_tor(end)=1;

  ne_si=spline(frp(:,1),frp(:,2),rho_pol)*transform_1013overcm3_to_oneoverm3;


  frp=load([path_to_shot, '/Te_ida_', shot_designation,'_rhopol.dat']);
  Te_eV=spline(frp(:,1),frp(:,2),rho_pol)*transform_keV_to_eV;


  frp=load([path_to_shot, '/Ti_cez_', shot_designation,'_rhopol.dat']);
  rho_fit2=rho_pol.^2;

  % octave does not have a 'fit' function, only 'fsolve'. Thus this was replaced.
  %[fitobject,gof] = fit(frp(:,1).^2, frp(:,2) ,'poly6');
  %fit2 = fitobject.p1.*(rho_fit2.^6) + fitobject.p2.*(rho_fit2.^5) + fitobject.p3.*(rho_fit2.^4) + fitobject.p4.*(rho_fit2.^3) + fitobject.p5.*(rho_fit2.^2) + fitobject.p6.*rho_fit2+fitobject.p7;
  [fitobject, gof] = polyfit(frp(:,1).^2, frp(:,2), 6);
  fit2 = polyval(fitobject, rho_fit2);

  Ti_eV=fit2*transform_keV_to_eV;

  if do_plots
    figure
    plot(rho_pol,Ti_eV*transform_eV_to_keV,'r',frp(:,1),frp(:,3),'b--',frp(:,1),frp(:,4),'b--')
    xlabel('\rho_{pol}')
    ylabel('T_i  [keV]')
  end


  frp=load([path_to_shot, '/vrot_cez_', shot_designation,'_rhopol.dat']);
  frr=load([path_to_shot, '/vrot_cez_', shot_designation,'_R.dat']);

  frp(:,2)=frp(:,2)./frr(:,1);
  frp(:,3)=frp(:,3)./frr(:,1);
  frp(:,4)=frp(:,4)./frr(:,1);

  frp(:,2:4)=frp(:,2:4)*transform_kHz_to_Hz;

  [fitobject, gof] = polyfit(frp(:,1).^2, frp(:,2), 6);
  fit2 = polyval(fitobject, rho_fit2);

  vrot=fit2;

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
