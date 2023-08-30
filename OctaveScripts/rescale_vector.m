% Rescale (and maybe interpolate) function/vector.
%
% input:
% ------
% x_grid_old: vector, old grid.
% y_old_grid_old: vector, old values on old grid.
% z_old_grid_old: vector which should be rescaled, old values on old grid.
% x_grid_new: vector, new grid. Might be empty, in this case it is
%   assumed that old and new grid are the same.
% y_new_grid_new: vector, the new values on the new grid.
%
% output:
% -------
% z_new
function [z_new] = rescale_vector(x_grid_old, y_old_grid_old, z_old_grid_old, x_grid_new, y_new_grid_new)
  % Make sure we have all required values on the new grid.
  if ~isempty(x_grid_new)
    % interpolation necessary
    y_old_grid_new = spline(x_grid_old, y_old_grid_old, x_grid_new);
    z_old_new_grid = spline(x_grid_old, z_old_grid_old, x_grid_new);
  else
    % no difference
    y_old_grid_new = y_old_grid_old;
    z_old_new_grid = z_old_grid_old;
  end
  z_new = z_old_new_grid .* y_new_grid_new ./ y_old_grid_new;
end
