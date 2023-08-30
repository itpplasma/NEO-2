% Get derivative of volume as function of normalized toroidal flux
%
% This function calculates the derivative of the volume as a function of
% the normalized toroidal magnetic flux from a neo-2 multispecies output
% file structure.
%
% \attenation Still has to be checked if it is correct!
%
% input:
% ------
% h5file: structure that contains the fields of a neo-2 multispecies
%   output file. The data of these is used to calculate dV/ds.
%
% output:
% -------
% dVds: array with the values of the derivative.
%
% sideeffects:
% ------------
% none
%
% limitations:
% ------------
% Required field names are hardcoded.
function dVds = get_dVds_in_SI(h5file)

  VOLUME_CGS_TO_SI = 1.0e-6;

  avb2 = h5file.avbhat2 .* h5file.Bref .* h5file.Bref;
  boozer_psi_pr = h5file.psi_pr_hat .* h5file.Bref;

  dVds = (2*pi)^2*(h5file.aiota .* h5file.bcovar_tht + h5file.bcovar_phi)./avb2*boozer_psi_pr;

  % as neo-2 uses cgs units, this should be in cm^3/?.
  dVds = dVds*VOLUME_CGS_TO_SI;

end
