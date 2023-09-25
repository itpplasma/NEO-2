% Generate script for a run which sets the required input paramters
%
% Note that you can run a script without starting starting matlab first
% by using:
% matlab -nodesktop -nosplash -r "path('$NEO2PATH/OctaveScripts', path); script_name; exit(0);"
%
% Input:
%   hdf5FileName: name for hdf5 file that is produced. Defaults to
%     'multi_spec_Valentin.in'.
%   path_to_shot: path to the shoot where the data is. Defaults to
%     'SHOT35568/'.
%   data_source: a structure, is passed to function 'load_profile_data',
%     see there for details.
%   species_definition: three dimensional array. The size of the first
%     dimension determines the number of radial points to use. The
%     second dimension determines the number of species. The third
%     dimension has two entries, refering to the charge and mass of the
%     corresponding species at the respecitve radial location.
%     Might be hard to imagine why the mass should vary, but the charge
%     could simply vary due to the change in temperature.
%   isw_Vphi_loc:
%   species_tag_Vphi: index of species for which the rotation velocity is
%     given.
%   input_unit_type: switch to determine what type/unit (and thus what
%     conversions have to be done) the input has.
%   bounds: two array index determining bounds for flux surface label in
%     which to put the points.
%   switch_grid: integer value (optional), passed to function
%     'load_profile_data'. See documentation of that function for more
%     details.
%
% output:
% -------
% None
%
% sideeffects:
% ------------
% Creates file to which the profiles are written.
%
% limitations:
% ------------
% Function is so far not capable of handling more than one ion species,
%   i.e. two species total.
function generate_neo2_profile(hdf5FileName, path_to_shot, data_source, species_definition, isw_Vphi_loc, species_tag_Vphi, input_unit_type, bounds, switch_grid)
  if nargin() < 1 || isempty(hdf5FileName)
    hdf5FileName = 'multi_spec_Valentin.in';
  end
  if nargin() < 2 || isempty(path_to_shot)
    path_to_shot = 'SHOT35568/';
  end
  if nargin() < 3 || isempty(data_source)
    shot_designation = '35568_t2.6881'
    data_source.rhopoloidal.filename = ['ne_ida_', shot_designation,'_rhopol.dat'];
    data_source.rhopoloidal.column = 1;
    data_source.rhotoroidal.filename = ['ne_ida_', shot_designation,'_rhotor.dat'];
    data_source.rhotoroidal.column = 1;
    data_source.electron_density.filename = ['ne_ida_', shot_designation,'_rhopol.dat'];
    data_source.electron_density.column = 2;
    data_source.electron_temperature.filename = ['Te_ida_', shot_designation,'_rhopol.dat'];
    data_source.electron_temperature.column = 2;
    data_source.ion_temperature.filename = ['Ti_cez_', shot_designation,'_rhopol.dat'];
    data_source.ion_temperature.column = 2;
    data_source.rotation_velocity.filename = ['vrot_cez_', shot_designation,'_rhopol.dat'];
    data_source.rotation_velocity.column = 2;
    data_source.major_radius.filename = ['vrot_cez_', shot_designation,'_R.dat'];
    data_source.major_radius.column = 1;
  end
  if nargin() < 4 || isempty(species_definition)
    Zi=1; % ion charge number (deuterium)
    mi=3.3436e-24; % ion mass (deuterium)
    Ze=-1;% electron charge number
    me=9.1094e-28; % electron mass
    species_definition(1,1,1) = Ze;
    species_definition(1,1,2) = me;
    species_definition(1,2,1) = Zi;
    species_definition(1,2,2) = mi;
    species_definition = repmat(species_definition, gridpoints, 1, 1);
  end
  if nargin() < 5 || isempty(isw_Vphi_loc)
    isw_Vphi_loc = 0;
  end
  if nargin() < 6 || isempty(species_tag_Vphi)
    species_tag_Vphi = 2;
  end
  if nargin() < 7 || isempty(input_unit_type)
    input_unit_type = 1;
  end
  if nargin() < 8 || isempty(bounds)
    gridpoints = size(species_definition, 1);
  else
    gridpoints = [size(species_definition, 1), bounds(1), bounds(2)];
  end
  % \attention NO isemtpy!, parameter is just passed to load_profile_data
  %   thus the function is responsible for default value. Here we make
  %   just sure the parameter exists.
  if nargin() < 9
    switch_grid = [];
  end
  num_species = size(species_definition, 2);

  %% UNIT CONV DEF
  ELEMENTARY_CHARGE_SI = 1.60217662e-19;
  SPEED_OF_LIGHT_SI = 2.99792458e8;
  EV_TO_SI = ELEMENTARY_CHARGE_SI;
  DENSITY_SI_TO_CGS = 1e-6;
  ENERGY_TO_CGS = 1e7;
  EV_TO_CGS = EV_TO_SI * ENERGY_TO_CGS;
  CHARGE_SI_TO_CGS = 10*SPEED_OF_LIGHT_SI;% Conversion factor is not speed of light, but 10c_si.

  [rho_pol, rho_tor, ne_si, Ti_eV, Te_eV, vrot] = load_profile_data(path_to_shot, data_source, gridpoints, 0, input_unit_type, switch_grid);

  %% BOOZER S
  boozer_s = (rho_tor).^2;

  %% PROFILES
  ne = ne_si * DENSITY_SI_TO_CGS;
  ni = ne;
  Te = Te_eV * EV_TO_CGS;
  Ti = Ti_eV * EV_TO_CGS;

  %% PROFILE DERIVATIVES
  derivate = @ (s,t) [(s(3)-s(1))./(t(3)-t(1)), (s(3:end)-s(1:end-2))./(t(3:end)-t(1:end-2)), (s(end)-s(end-2))./(t(end)-t(end-2))];

  dTe_ov_ds = derivate(Te,boozer_s);
  dTi_ov_ds = derivate(Ti,boozer_s);
  dne_ov_ds = derivate(ne,boozer_s);
  dni_ov_ds = derivate(ni,boozer_s);

  %% SPEC INIT
  num_radial_pts = length(rho_pol);
  species_tag=1:num_species;

  rel_stages = zeros(1, num_radial_pts);
  rel_stages(:) = species_tag_Vphi;

  %% KAPPA
  % Coulomb logarithm (set as species-independent - see E A Belli and J Candy PPCF 54 015015 (2012))
  log_Lambda = 39.1 - 1.15*log10(ne_si) + 2.3*log10(Te_eV*1e-3);

  % Conversion cgs-units
  e_e = -ELEMENTARY_CHARGE_SI * CHARGE_SI_TO_CGS;%[statC]
  e_i = +ELEMENTARY_CHARGE_SI * CHARGE_SI_TO_CGS;%[statC]
  Te_ov_ee = -Te/e_e;%[erg/statC]
  Ti_ov_ei = Ti/e_i;%[erg/statC]

  % mean free path
  lc_e=(3/(4*sqrt(pi)))*(Te_ov_ee.^2) ./ (ne.*(e_e^2).*log_Lambda);
  lc_i=(3/(4*sqrt(pi)))*(Ti_ov_ei.^2) ./ (ni.*(e_i^2).*log_Lambda);

  % Compute kappa for electrons and main ion species
  kappa_e = 2./lc_e;
  kappa_i = 2./lc_i;


  %% WRITE HDF5
  % Explicitly remove existing file, to avoid problems when writing.
  system(['rm ',hdf5FileName]);

  % The part below does not work in octave, as h5create is missing (up to at least 7.3.0).
  h5create(hdf5FileName,'/num_radial_pts',1,'DataType','int32')
  h5write(hdf5FileName,'/num_radial_pts',int32(num_radial_pts))
  h5create(hdf5FileName,'/num_species',1,'DataType','int32')
  h5write(hdf5FileName,'/num_species',int32(num_species))

  h5create(hdf5FileName,'/species_tag',num_species,'DataType','int32')
  h5write(hdf5FileName,'/species_tag',int32(species_tag))
  h5create(hdf5FileName,'/species_def',[num_radial_pts,num_species, 2],'DataType','double')
  h5write(hdf5FileName,'/species_def',species_definition)
  h5writeatt(hdf5FileName, '/species_def', 'unit', 'e; g');


  h5create(hdf5FileName,'/boozer_s',num_radial_pts,'DataType','double')
  h5write(hdf5FileName,'/boozer_s',boozer_s)
  h5writeatt(hdf5FileName, '/boozer_s', 'unit', '1');

  h5create(hdf5FileName,'/rho_pol',num_radial_pts,'DataType','double')
  h5write(hdf5FileName,'/rho_pol',rho_pol)
  h5writeatt(hdf5FileName, '/rho_pol', 'unit', '1');

  h5create(hdf5FileName,'/Vphi',num_radial_pts,'DataType','double')
  h5write(hdf5FileName,'/Vphi',vrot)
  h5writeatt(hdf5FileName, '/Vphi', 'unit', 'rad / s');
  h5create(hdf5FileName,'/species_tag_Vphi',1,'DataType','int32')
  h5write(hdf5FileName,'/species_tag_Vphi',int32(species_tag_Vphi))
  h5create(hdf5FileName,'/isw_Vphi_loc',1,'DataType','int32')
  h5write(hdf5FileName,'/isw_Vphi_loc',int32(isw_Vphi_loc))

  h5create(hdf5FileName,'/rel_stages',num_radial_pts,'DataType','int32')
  h5write(hdf5FileName,'/rel_stages',int32(rel_stages))

  h5create(hdf5FileName,'/T_prof',[num_radial_pts,num_species],'DataType','double')
  h5write(hdf5FileName,'/T_prof',[Te;Ti]')
  h5writeatt(hdf5FileName, '/T_prof', 'unit', 'erg');
  h5create(hdf5FileName,'/dT_ov_ds_prof',[num_radial_pts,num_species],'DataType','double')
  h5write(hdf5FileName,'/dT_ov_ds_prof',[dTe_ov_ds; dTi_ov_ds]')

  h5create(hdf5FileName,'/n_prof',[num_radial_pts,num_species],'DataType','double')
  h5write(hdf5FileName,'/n_prof',[ne; ni]')
  h5writeatt(hdf5FileName, '/n_prof', 'unit', '1 / cm^3');
  h5create(hdf5FileName,'/dn_ov_ds_prof',[num_radial_pts,num_species],'DataType','double')
  h5write(hdf5FileName,'/dn_ov_ds_prof',[dne_ov_ds; dni_ov_ds]')

  h5create(hdf5FileName,'/kappa_prof',[num_radial_pts,num_species],'DataType','double')
  h5write(hdf5FileName,'/kappa_prof',[kappa_e; kappa_i]')
end
