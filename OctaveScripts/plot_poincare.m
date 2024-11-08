% Make a poincare plot of data from magnetic.h5
%
% input:
% ------
% filename: string, name (and path) of file to load.
%
% output:
% -------
% x1,x2,x3: arrays, 'slices' of the corresponding fields of 'fieldperiod'
%   in the magnetic.h5 file (1 and 3 are those actually used for the
%   plot).
%
% side effects:
% -------------
% Makes a plot
function [x1, x2, x3] = plot_poincare(filename)
  %~ h5file = '/path/to/Neo2/RunsByDate/2017_07_Bootstrap/w7as/conl_over_mfp=1m1-boozer_s=0.25d0-mag_nperiod_min=100-bsfunc_local_err=3m2/magnetics.h5';
  %~ h5 = h52struct(h5file);
  h5 = load(filename);

  fnames = fieldnames(h5.fieldperiod);
  idx = 1;
  for n = 1:numel(fnames)
    x1(n) = h5.fieldperiod.(fnames{n}).x1(idx);
    x2(n) = h5.fieldperiod.(fnames{n}).x2(idx);
    x3(n) = h5.fieldperiod.(fnames{n}).x3(idx);
  end

  plot(x1,x3, '.')
end
