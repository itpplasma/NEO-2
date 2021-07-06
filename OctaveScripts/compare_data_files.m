#! /usr/bin/octave -q

% script for comparing two data files.
%
% expects four arguments, names(+path) of the two files, absolute
% accuracy, and relative accuracy.
%
% Unfortunately there seems to be no way to export the result.

if nargin < 4
  error("not enough input arguments!")
end

arg_list = argv();
file_1 = arg_list{1};
file_2 = arg_list{2};
abs_accuracy = str2double(arg_list{3});
rel_accuracy = str2double(arg_list{4});

try
  load(file_1);
catch
  error(['first file (', file_1, ') could not be loaded']);
end

try
  load(file_2);
catch
  error(['second file (', file_1, ') could not be loaded']);
end

if (ndims(file_1) != ndims(file_2)) || (size(file_1) != size(file_2))
  diff = abs(file_1 -file_2);

  m = max(file_1(:));

  r = all(diff < abs_accuracy) & all(diff ./ m < rel_accuracy);
end
