% Caclulate ripple plateau coefficient in cgs units for given neo2 file.
%
% Reads a neo-2 output file, given by name(+path) und uses the necessary
% quantities to calculate the ripple plateau coefficient in cgs units
% for given toroidal mode number and species. Mode number enters as
% factor, while species enters via mass and temperature.
%
% input:
% ------
% filename: string, name+path of the file for which to do the calculations.
% n: integer, toroidal mode number.
% species: integer, for which species to calculate the coefficient.
%
% output:
% -------
% Drp: (row/column?) vector, with same number of elements as for radial
%   direction. The ripple plateau coefficient in cgs units as function
%   of radial direction.
function [Drp] = calculate_cgs_ripple_plateau_coefficient_for_file(filename, n, species)
  statC_to_coulomb = 3.33564e-10;
  ms_to_cms = 1.0e+2; % meant is m/s and cm/s

  c = 2.99792458e+8*ms_to_cms; % cm/s
  e0 = 1.602176634e-19 ./ statC_to_coulomb; % statC


  h = load('neo2_multispecies_out.h5');

  F = 1.0 ./ h.avnabpsi.^2 ./ (h.psi_pr_hat.^2 .*h.Bref) .*(h.aiota.*h.bcovar_tht + h.bcovar_phi) ./ h.aiota.^2;

  omega_co = e0*h.Bref ./ (h.m_spec(species,:)*c);
  vthermal = sqrt(2*h.T_spec(species,:) ./ (h.m_spec(species,:)));
  rho_L = vthermal ./ omega_co;

  Drp = sqrt(pi)/4*n *vthermal.*(rho_L.^2) .*F.*h.eps_M_2;
end
