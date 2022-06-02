% Convert in- and output of neo-2 to neo-rt input.
%
% input:
% ------
% infileneo2out: string, name of neo-2 output file to read.
% outfilename: string, name under which to save the results.
% number_surfaces: number of flux surfaces to evaluate/write to file.
%
% output:
% -------
% none
%
% sideeffects:
% ------------
% Creates/overwrites file.
function neo2_to_neort(infileneo2out, outfilename, number_surfaces)
    nout = load(infileneo2out);

  ds = 1.0 ./ number_surfaces;

  f = fopen(outfilename, 'w');
  for k = 1:number_surfaces
    s = ds*k;

    mt = spline(nout.boozer_s, nout.MtOvR(2,:), s);
    n = spline(nout.boozer_s, nout.n_spec(2,:), s);
    fprintf(f, '%13.7f  %13.7f  %13.7f\n', s, mt, n);
  end
  fclose(f);
end
