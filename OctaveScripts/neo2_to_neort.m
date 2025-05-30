% Convert in- and output of neo-2 to neo-rt input.
%
% input:
% ------
% infileneo2out: string, name of neo-2 output file to read.
% outfilename_prefix: string, prefix for name under which to save the
%     results. E.g. if outfilename_prefix='test_' the output files will
%     be 'test_profile.in' and test_plasma.in'.
% number_surfaces: number of flux surfaces to evaluate/write to file.
%
% output:
% -------
% none
%
% sideeffects:
% ------------
% Creates/overwrites files.
function neo2_to_neort(infileneo2out, outfilename_prefix, number_surfaces)
  nout = load(infileneo2out);

  ds = 1.0 ./ number_surfaces;

  cm3_to_m3 = 1.0e+6;
  ev_to_cgs = 1.6022e-12;
  kg_to_dalton = 1.0/1.66053906660e-27;
  g_to_kg = 1.0e-3;

  species = 2;

  s = ds*[1:number_surfaces];
  mt = spline(nout.boozer_s, nout.MtOvR(species,:) .* nout.R0, s);
  vth = spline(nout.boozer_s, sqrt(2*nout.T_spec(species,:) ./ nout.m_spec(species, :)), s);

  f = fopen([outfilename_prefix, 'profile.in'], 'w');
  for k = 1:number_surfaces
    fprintf(f, '%13.7f  %13.7f  %13.7f\n', s(k), mt(k), vth(k));
  end
  fclose(f);

  n1 = spline(nout.boozer_s, nout.n_spec(species,:), s);
  n2 = zeros(1, number_surfaces);

  T1 = spline(nout.boozer_s, nout.T_spec(species,:)/ev_to_cgs, s);
  T2 = T1;
  Te = spline(nout.boozer_s, nout.T_spec(1,:)/ev_to_cgs, s);

  f = fopen([outfilename_prefix, 'plasma.in'], 'w');
  fprintf(f, '%% N am1 am2 Z1 Z2\n');
  fprintf(f, '%4i %9.3f %9.3f %9.3f %9.3f\n',
              number_surfaces,
              nout.m_spec(species,1)*g_to_kg*kg_to_dalton,
              nout.m_spec(species,1)*g_to_kg*kg_to_dalton,
              nout.z_spec(species,1),
              nout.z_spec(species,1));
  fprintf(f, '%% s ni_1[cm^-3] ni_2[cm^-3] Ti_1[eV] Ti_2[eV] Te[eV]\n');
  for k = 1:number_surfaces
    fprintf(f, '%13.7e  %13.7e  %13.7e  %13.7e  %13.7e  %13.7e\n', ...
        s(k), n1(k), n2(k), T1(k), T2(k), Te(k));
  end
  fclose(f);
end
