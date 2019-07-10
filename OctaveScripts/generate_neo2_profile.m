% matlab -nodesktop -nosplash -r "path('/afs/itp.tugraz.at/user/buchholz/Programs/neo-2/OctaveScripts', path); generate_neo2_profile; exit(0);"

clear all %#ok<CLALL>

%% INPUT

hdf5FileName = 'multi_spec_Valentin.in';

gridpoints = 100;
num_species = 2;
isw_Vphi_loc = 0;
species_tag_Vphi = 2;

[rho_pol, rho_tor, ne_si, Ti_eV, Te_eV, vrot] = load_profile_data('SHOT35568/', '35568_t2.6881', gridpoints, 0);

%% UNIT CONV DEF
ELEMENTARY_CHARGE = 1.60217662e-19;
SPEED_OF_LIGHT = 2.99792458e9;
EV_TO_SI = ELEMENTARY_CHARGE;
DENSITY_TO_CGS = 1e-6;
SPEED_TO_CGS = 100;
ENERGY_TO_CGS = 1e7;
VOLT_TO_CGS = ENERGY_TO_CGS/SPEED_OF_LIGHT;
EV_TO_CGS = EV_TO_SI * ENERGY_TO_CGS;
CHARGE_TO_CGS = SPEED_OF_LIGHT;

%% SPEC DEF
Zi=1; % ion charge number (deuterium)
mi=3.3436e-24; % ion mass (deuterium)
Ze=-1;% electron charge number
me=9.1094e-28; % electron mass

%% BOOZER S
boozer_s = (rho_tor).^2;

%% PROFILES
ne = ne_si * DENSITY_TO_CGS;
ni = ne;
Te = Te_eV * EV_TO_CGS;
Ti = Ti_eV * EV_TO_CGS;

%% PROFILE DERIVATIVES
derivate = @ (s,t) [(s(3)-s(1)/t(3)-t(1)), (s(3:end)-s(1:end-2))./(t(3:end)-t(1:end-2)), (s(end)-s(end-2)/t(end)-t(end-2))];

dTe_ov_ds = derivate(Te,boozer_s);
dTi_ov_ds = derivate(Ti,boozer_s);
dne_ov_ds = derivate(ne,boozer_s);
dni_ov_ds = derivate(ni,boozer_s);

%% SPEC INIT
num_radial_pts = length(rho_pol);
species_tag=1:num_species;
species_def=zeros(num_radial_pts,num_species,2);

species_def(:,1,1) = Ze;
species_def(:,1,2) = me;
species_def(:,2,1) = Zi;
species_def(:,2,2) = mi;

rel_stages = zeros(1, num_radial_pts);
rel_stages(:) = species_tag_Vphi;

%% KAPPA
% Coulomb logarithm (set as species-independent - see E A Belli an J Candy PPCF 54 015015 (2012))
log_Lambda = 39.1 - 1.15*log10(ne_si) + 2.3*log10(Te_eV*1e-3);

% Conversion cgs-units
e_e = -ELEMENTARY_CHARGE * CHARGE_TO_CGS;%[statC]
e_i = +ELEMENTARY_CHARGE * CHARGE_TO_CGS;%[statC]
Te_ov_ee = -Te/e_e;%[erg/statC]
Ti_ov_ei = Ti/e_i;%[erg/statC]

% mean free path
lc_e=(3/(4*sqrt(pi)))*(Te_ov_ee.^2) ./ (ne.*(e_e^2).*log_Lambda);
lc_i=(3/(4*sqrt(pi)))*(Ti_ov_ei.^2) ./ (ni.*(e_i^2).*log_Lambda);

% Compute kappa for electrons and main ion species
kappa_e = 2./lc_e;
kappa_i = 2./lc_i;


%% WRITE HDF5
system(['rm ',hdf5FileName])

% The part below does not work in octave (at least 4.0.3), as h5create is missing.
h5create(hdf5FileName,'/num_radial_pts',1,'DataType','int32')
h5write(hdf5FileName,'/num_radial_pts',int32(num_radial_pts))
h5create(hdf5FileName,'/num_species',1,'DataType','int32')
h5write(hdf5FileName,'/num_species',int32(num_species))

h5create(hdf5FileName,'/species_tag',num_species,'DataType','int32')
h5write(hdf5FileName,'/species_tag',int32(species_tag))
h5create(hdf5FileName,'/species_def',[num_radial_pts,num_species, 2],'DataType','double')
h5write(hdf5FileName,'/species_def',species_def)


h5create(hdf5FileName,'/boozer_s',num_radial_pts,'DataType','double')
h5write(hdf5FileName,'/boozer_s',boozer_s)

h5create(hdf5FileName,'/rho_pol',num_radial_pts,'DataType','double')
h5write(hdf5FileName,'/rho_pol',rho_pol)

h5create(hdf5FileName,'/Vphi',num_radial_pts,'DataType','double')
h5write(hdf5FileName,'/Vphi',vrot)
h5create(hdf5FileName,'/species_tag_Vphi',1,'DataType','int32')
h5write(hdf5FileName,'/species_tag_Vphi',int32(species_tag_Vphi))
h5create(hdf5FileName,'/isw_Vphi_loc',1,'DataType','int32')
h5write(hdf5FileName,'/isw_Vphi_loc',int32(isw_Vphi_loc))

h5create(hdf5FileName,'/rel_stages',num_radial_pts,'DataType','int32')
h5write(hdf5FileName,'/rel_stages',int32(rel_stages))

h5create(hdf5FileName,'/T_prof',[num_radial_pts,num_species],'DataType','double')
h5write(hdf5FileName,'/T_prof',[Te;Ti]')
h5create(hdf5FileName,'/dT_ov_ds_prof',[num_radial_pts,num_species],'DataType','double')
h5write(hdf5FileName,'/dT_ov_ds_prof',[dTe_ov_ds; dTi_ov_ds]')

h5create(hdf5FileName,'/n_prof',[num_radial_pts,num_species],'DataType','double')
h5write(hdf5FileName,'/n_prof',[ne; ni]')
h5create(hdf5FileName,'/dn_ov_ds_prof',[num_radial_pts,num_species],'DataType','double')
h5write(hdf5FileName,'/dn_ov_ds_prof',[dne_ov_ds; dni_ov_ds]')

h5create(hdf5FileName,'/kappa_prof',[num_radial_pts,num_species],'DataType','double')
h5write(hdf5FileName,'/kappa_prof',[kappa_e; kappa_i]')
