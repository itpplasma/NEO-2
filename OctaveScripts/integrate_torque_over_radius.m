% function [integrated_torque, s_grid] = integrate_torque_over_radius(h5file, exclude_region)
%
% Return integrated torque as function of radius, together with grid
% over which was integrated and logical array to recreate the grid from
% the original grid.
%
% input:
% ------
% h5file: structure, as returned from a 'load(h5file)'.
% exclude_region: optional 2 element vector [lower_bound, upper_bound].
%   If present determines region to exclude from the computation of
%   torque. Only values of s either smaller as the first value, or
%   larger as the second value are used. Default value ist [2, 2].
% ignore_large_values: if true than all points for which the magnitude
%   of electron torque is larger than 1.0e+10 are ignored (in addition
%   to the exclude region). If not given or false, then no further
%   restricions are made. Defaults to false.
%
% output:
% -------
% integrated_torque: vector with the torque values.
% s_grid: s_grid of torque values, might differ from input due to
%   exclude_region.
% L: logical array, which determines the radial points over which is
%   integrated.
function [integrated_torque, s_grid, L] = integrate_torque_over_radius(h5file, exclude_region, ignore_large_values)
  if nargin < 2 || isempty(exclude_region)
    exclude_region = [2, 2];
  end

  if nargin < 3 || isempty(ignore_large_values)
    ignore_large_values = false;
  end

  torque_cgs_to_si = 1.0e-7;

  L1 = h5file.boozer_s < exclude_region(1);
  L2 = h5file.boozer_s > exclude_region(2);

  % Handling of large values in torque.
  L3 = true(size(h5file.boozer_s));
  if (ignore_large_values)
    L3 = abs(h5file.TphiNA_spec(1, :)) < 1.0e+10;
  end
  L1 = L1 & L3;
  L2 = L2 & L3;
  L = L1 | L2;

  s_grid1 = h5file.boozer_s(L1);
  s_grid2 = h5file.boozer_s(L2);
  s_grid = [s_grid1, s_grid2];

  avb2 = h5file.avbhat2 .* h5file.Bref .* h5file.Bref;
  boozer_psi_pr = h5file.psi_pr_hat .* h5file.Bref;

  for k = 1:size(h5file.TphiNA_spec, 1)
    if any(L1)
      integrated_torque_temp = torque_cgs_to_si*(4*pi*pi)*boozer_psi_pr(end) ...
        .* cumtrapz(s_grid1,...
          h5file.TphiNA_spec(k, L1).*(h5file.aiota(L1) .* h5file.bcovar_tht(L1) + h5file.bcovar_phi(L1))./avb2(L1));
      offset = integrated_torque_temp(end);
    else
      offset = 0.0;
      integrated_torque_temp = [];
    end

    if any(L2)
      integrated_torque(k,:) = [integrated_torque_temp,...
        torque_cgs_to_si*(4*pi*pi)*boozer_psi_pr(end) ...
        .* cumtrapz(s_grid2,...
          h5file.TphiNA_spec(k,L2).*(h5file.aiota(L2) .* h5file.bcovar_tht(L2) + h5file.bcovar_phi(L2))./avb2(L2))...
         + offset];
    else
      integrated_torque(k,:) = integrated_torque_temp;
    end
  end
  %TphiNA_int_io = (4*pi*pi)*boozer_psi_pr[-1]*integrate.cumtrapz(array(TphiNA_io)*(array(aiota)*array(bcovar_tht)+array(bcovar_phi))/avb2, boozer_s, 0)
end
