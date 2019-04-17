% function [output_data] = create_displacement(filename, nr_points_rho, rho_q1, r_res, ksi22_max, tearing_parameter, plot_data)
%
% \brief Displacement function for toroidal perturbation modes.
%
% Default values correspond to the value used for 35568, t=2.68.
%
% Input:
% ------
% filename: name of file where to store the calculated values. [displacement_morepoints.txt]
% nr_points_rho: number of points to use for the radial, equidistant grid. [1000]
% rho_q1: SXR fishbones [0.25]
% r_res: resonant surface for (3,2), i.e. q(r_res)=3/2. [0.44945]
% ksi22_max: amplitude for the modes. [5]
% tearing_parameter: parameters of the function used to calculate the
%   displacement for the tearing mode. [0.06, 1, 1, 1]
% plot_data: logical, if true the calculated displacements are plotted
%   against rho. [false]
%
% Output:
% -------
% output_data: 2d array, where the columns are rho, rho^2, displacement
%   for the kink and tearing mode, respectively. Number of rows is equal
%   to nr_points_rho.
function [output_data] = create_displacement(filename, nr_points_rho, rho_q1, r_res, ksi22_max, tearing_parameter, plot_data)
  % define main parameters
  %~ rho_q1 = 0.25; %SXR fishbones rho_tor = 0.2/0.24 rho_pol=0.33-0.29
  %~ r_res = 0.44945; % resonant surface for (3,2) rho_tor
  %~ ksi22_max = 5;
  %~ W = 0.012 % width of the island in rho_tor
  if nargin < 1 || isempty(filename)
    filename = 'displacement_morepoints.txt';
  end
  if nargin < 2 || isempty(nr_points_rho)
    nr_points_rho = 1000;
  end
  if nargin < 3 || isempty(rho_q1)
    rho_q1 = 0.25;
  end
  if nargin < 4 || isempty(r_res)
    r_res = 0.44945;
  end
  if nargin < 5 || isempty(ksi22_max)
    ksi22_max = 5;
  end
  if nargin < 6 || isempty(tearing_parameter)
    tearing_parameter = [0.06, 1, 1, 1];
  end
  if nargin < 7 || isempty(plot_data)
    plot_data = false;
  end

  % Define radial grid.
  rho = linspace(0, 1, nr_points_rho);
  rho2 = rho.^2;

  % Construct function of the displacement for (2,2)
  displacement_kink = create_displacement_kink(rho, ksi22_max, rho_q1);

  % Construct function of the displacement for (3,2)
  displacement_tearing = create_displacement_tearing(rho, r_res, ksi22_max, tearing_parameter);

  output_data = ([rho(:), rho2(:), displacement_kink(:), displacement_tearing(:)]);

  save(filename, 'output_data', '-ascii');

  if plot_data
    figure(1)
    plot(rho, displacement_kink, rho, displacement_tearing)
    grid;
    legend('(2,2)','(3,2)')
    xlabel('\rho_{tor}')
  end
end

function [displacement_kink] = create_displacement_kink(rho, ksi22_max, rho_q1)
  a = -4*ksi22_max/(rho_q1^2);
  b = 4*ksi22_max/rho_q1;

  displacement_kink = a.*(rho.^2) + b.*rho;
  displacement_kink(displacement_kink<0) = 0;
end

function [displacement_tearing] = create_displacement_tearing(rho, rho_res, amplitude, tearing_parameter)
  delta_rho = rho-rho_res;

  ksi32inner = -(delta_rho.^tearing_parameter(2)).*exp(-(abs(delta_rho)./tearing_parameter(1)).^tearing_parameter(3));
  ksi32outer = -(delta_rho.^tearing_parameter(2)).*exp(-(abs(delta_rho)./tearing_parameter(1)).^tearing_parameter(3))*tearing_parameter(4);

  ksi32inner(ksi32inner > 0) = 0;
  ksi32outer(ksi32outer <= 0) = 0;

  displacement_tearing = ksi32inner + ksi32outer;
  maxdisplacement_tearing = max(displacement_tearing);

  % Renormalize, to set the amplitude.
  displacement_tearing = amplitude*displacement_tearing/maxdisplacement_tearing;
end
