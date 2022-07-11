% Convert zero-based index to a folder name.
%
% Intendet for example to convert the index of a condor run (e.g. from
% an error file like err.158), to a corresponding folder where the
% actual run was done.
%
% Example usage:
%
% folder = get_name_from_zero_based_index([158, 171; 172, 177], 'es_*')
% folder =
% {
%   [1,1] = es_0p05532
%   [2,1] = es_0p06591
%   [3,1] = es_0p06512
%   [4,1] = es_0p06996
% }
%
% input:
% ------
% index: zero based index of the folder. May also be a vector or a
%   matrix to get the folder name for multiple values at once.
% folders: string which determines the files for which the index is.
%   Usually this should be a regular expression.
%   Defaults to 'es_*', i.e. all files which start with 'es_' in the
%   current folder.
%
% output:
% -------
% cell array which contains the name of the folder(s) belonging to the
%   given index(indices).
%
% sideeffects:
% ------------
% none
%
% limitations:
% ------------
% Works only for folders, not for files as '-d' is used as option to ls.
function folder = get_name_from_zero_based_index(index, folders)
  if nargin < 2 || isempty(folders)
    folders = 'es_*';
  end

  % Create the list of folders.
  d = ls(['-d ', folders]);
  % Select those with matching index (+1 because octave uses one-based
  % indices).
  folder = mat2cell(d(index(:)+1, :), ones(size(index(:), 1),1));
end
