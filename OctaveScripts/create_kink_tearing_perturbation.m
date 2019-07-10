% Function for the creation of displacement and perturabtion files.
%
% This function will create a displacement file and three boozer files.
% The boozer files will have kink mode, tearing mode and both.
%
% input:
% ------
% file_base: name of the file (without ending, assumed to be '.bc'),
%   from which to read the equilibrium part of the field. Also used for
%   the name of the output file with perturbation.
% filename_displacement: Name under which the file with the
%   displacements is stored (including file ending).
% filename_add: string that might be added to the filename of the files
% for only kink/tearing modes.
% number_points: the number of points to use for the radial direction.
% function_parameters: array with seven values, used for the functions
%   that define the kink/tearing offset. This does not include amplitude
%   and phase.
% amplitudes: the amplitudes to use for the modes in the perturbation
%  files.
% phases: the phases to use for the modes in the perturbation files.
% plot_data: logical, if true a plot of the results is made.
function create_kink_tearing_perturbation(file_base, filename_displacement, filename_add,...
  number_points, function_parameters, amplitudes, phases, plot_data)

  create_displacement(filename_displacement,...
    number_points, function_parameters(1), function_parameters(2), function_parameters(3), function_parameters(4:7), plot_data);

  create_asdex_perturb(file_base, filename_displacement, amplitudes, phases, plot_data);

  amplitudes_ = [amplitudes(1), 0];
  create_asdex_perturb(file_base, filename_displacement, amplitudes_, phases,...
    plot_data, [file_base, '-pert_kink-n2-phase_kt-1pi', filename_add, '.bc']);

  amplitudes_ = [0, amplitudes(2)];
  create_asdex_perturb(file_base, filename_displacement, amplitudes_, phases,...
    plot_data, [file_base, '-pert_tear-n2-phase_kt-1pi', filename_add, '.bc']);
end
