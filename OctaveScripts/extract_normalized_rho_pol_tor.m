% Extract normalized poloidal and toroidal flux from efit_to_boozer.x output
%
% Extract normalized poloidal flux and normalized toroidal flux from
% the (additional) output of efit_to_boozer.x.
%
% Might be needed for the generation of profiles. The matlab function
% requires either one quantity as function of both or a file with booth
% quantities.
%
% Example usage:
%   Using default values for filenames.
%   [rho_pol, rho_tor] = extract_normalized_rho_pol_tor()
%
%   Not writing output file, note that string for outputfilename is not
%   empty. Using default value for inputfilename.
%   [rho_pol, rho_tor] = extract_normalized_rho_pol_tor([], ' ')
%
% input:
% ------
% inputfilename: name of the file from which to read the data.
%   Defaults to 'flux_functions.dat'.
% outputfilename: name of the file in which to store the normalized
%   data.
%   Note that the string is trimmed before being used to open the file.
%   If the string then is empty, no file will be created.
%   Defaults to 'rho_tor_vs_rho_pol.dat'.
%
% output:
% -------
% rho_pol: column vector, with normalized rho_pol values.
% rho_tor: column vector, with normalized rho_tor values.
%
% sideeffects:
% ------------
% Creates file if outputfilename does not contain only whitespace.
function [rho_pol, rho_tor] = extract_normalized_rho_pol_tor(inputfilename, outputfilename)
  if nargin < 1 || isempty(inputfilename)
    inputfilename = 'flux_functions.dat';
  end
  if nargin < 2 || isempty(outputfilename)
    outputfilename = 'rho_tor_vs_rho_pol.dat';
  end

  dat = load(inputfilename);

  rho_pol = dat(:, 4);
  rho_pol = rho_pol ./ rho_pol(end);

  rho_tor = dat(:, 6);
  rho_tor = rho_tor ./ rho_tor(end);

  if not(strcmp(strtrim(outputfilename), ''))
    fid = fopen(outputfilename, 'w');
    for k = 1:numel(rho_tor)
      fprintf(fid,' %11.4E %11.4E\n',[rho_pol(k), rho_tor(k)]);
    end
    fclose(fid);
  end

end
