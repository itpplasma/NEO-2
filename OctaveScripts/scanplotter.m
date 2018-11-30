#!/usr/bin/env octave

% Example script to show how one can read the data from a neo-2, using
% octave (based on a matlab version).
% Defines some functions, that might be of more general interest.
%
% Note that the h5read function should not be used with large hdf5
% files, as it will read in the whole file, and then select the
% requested part, i.e. for large files much memory is needed. Not sure
% how the size scales with representation in memory, but I would guess
% at least one to one, with tendency to larger memory representations.
%
% Usage of this script requires to set up the scan folders with specific
% name structure.

%collision parameter scan for tokamak
project_out='/temp/aradi_m/tokamak_collpar_scan_full';

file_name = 'fulltransp.h5';
file_name2 = 'neo2_config.h5';

param_value_array1={'1m5' '2m5' '5m5' '1m4' '2m4' '5m4' '1m3' '2m3' '5m3' '1m2' '2m2' '5m2' '1m1' '2m1' '5m1' '1' '2' '5' '10' '20' '50' '100' '200' '500' '1000'};
param_value_array2={'0p01' '0p09' '0p25' '0p49' '0p81' '0p99'};

% Get a field/part from a structure.
% It requires the structure itself and the 'path' to the desired
% field/part, where levels are separated by '/' or '.'.
% Consider for example the structure (intendation depicts level, i.e.
% n.a)
%
%   n
%     a
%       b
%       c
%     d
%     e
%     z
%       x
%       y
% here b,c,d,e,x and y may be arbitrary values (e.g. integer, float or
% array).
% Now the result of
% get_field_with_path_from_struct(n, '')
%   would be n
% get_field_with_path_from_struct(n, 'a/b')
%   would be b
% get_field_with_path_from_struct(n, 'a')
%   would be the structure a
function [f] = get_field_with_path_from_struct(struct, path)
  % For simplicity remove leading and trailing blanks.
  trimpath = strtrim(path);
  % Split the path+name of the object into single parts.
  % Repeated occurences of '/' should be ignored/removed, i.e. they do
  % not lead to an empty entry, according to a short test and according
  % to examples in documentation.
  temp = struct;
  cell_list = strsplit(trimpath, {'/', '.'});
  for i = 1:length(cell_list)
    % There should be no empty entries inbetween, as multiple delimeters
    % are merged by default, but this is not the case at the beginning
    % and the end, i.e. if the path starts/ends with a delimiter.
    if (length(cell_list{i}) > 0)
      tempstr = strjoin({'temp = temp.',cell_list{i} ,';'},'');
      eval(tempstr);
    end
  end

  f = temp;
end

% Function overloading with different number of arguments does not seem
% to be supported (type might be possible, at least with v4.4+).
%~ % As above, but with the fieldname as seperate parameter. This will just
%~ % Joind path and fieldname and pass the result to the version above.
%~ function [f] = get_field_with_path_from_struct(struct, path, fieldname)
  %~ f = get_field_with_path_from_struct(struct, strjoin({path, fieldname}, '.'));
%~ end

% Circumvent the point that octave (at least in version 4.0.3) has no
% h5read function, by defining our own.
% \Attention This function will read in the whole file, and then select
%   just the appropriate part.
%   Thus is is not suitable for large files.
function [h] = h5read(pathandfilename, pathandvariablename)
  t = load('-hdf5', pathandfilename);

  h = get_field_with_path_from_struct(t, pathandvariablename);
end

param_name1='conl_over_mfp';
param_path1='';
param_name2='boozer_s'
param_path2='settings';
param_elements1=length(param_value_array1);
param_elements2=length(param_value_array2);

for i1=1:param_elements1
  for i2=1:param_elements2
    folder_name = [project_out, '/', param_name1, '=', param_value_array1{i1}, '&', param_name2, '=', param_value_array2{i2}];
    disp(['Reading ', folder_name,'/',file_name])

    param_set1 = h5read([folder_name,'/',file_name], [param_path1, '/', param_name1]);

    param_set2 = h5read([folder_name,'/',file_name2], [param_path2, '/', param_name2]);

    gamma_out = h5read([folder_name,'/',file_name], 'gamma_out');

    conl_over_mfp(i1,i2) = param_set1;
    boozer_s(i1,i2) = param_set2;
    radial_transport_coeff(i1,i2) = gamma_out(1,1);
  end
end

loglog(conl_over_mfp, radial_transport_coeff, 'x-')
legend(param_value_array2)
