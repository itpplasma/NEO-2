% Scipt to plot bootstrap current.

%% Preamble
clear
addpath('/proj/plasma/Neo2/Interface/Matlab/')

%% Define filters
%filter_param = {};
%filter_value = {};
filter_param = {'leg', 'lag'};
filter_value = {1, 1};


%% Filter dirs
dirs = dir('.');
filtered_dirs_num = 0;
dirs_filtered = cell(0);
for k = 1:numel(dirs)
  if dirs(k).isdir
    dirname = dirs(k).name;

    if (exist([dirname, '/neo2_config.h5'], 'file') == 2) && (exist([dirname, '/efinal.h5'], 'file') == 2) ...
          && (exist([dirname, '/NOT_READY'], 'file') ~= 2)
      neo2config = h52struct([dirname, '/neo2_config.h5']);
      neo2configfields = fieldnames(neo2config);

      filter_ok = true;
      for l = 1:numel(neo2configfields)
        for f = 1:numel(filter_param)
          try
            if neo2config.(neo2configfields{l}).(filter_param{f}) ~= filter_value{f}
              filter_ok = false;
              break
            end
          end
        end
      end
      if (filter_ok)
        filtered_dirs_num = filtered_dirs_num + 1;
        dirs_filtered{filtered_dirs_num} = dirname;
      end

    end

  end
end

%% Plotter
for k = 1:numel(dirs_filtered)
  dirname = dirs_filtered{k};
  efinal      = h52struct([dirname, '/efinal.h5']);
  fulltransp  = h52struct([dirname, '/fulltransp.h5']);
  neo2config  = h52struct([dirname, '/neo2_config.h5']);

  mag_magfield(k)  = neo2config.settings.mag_magfield;
  mag_magfield_multiplier = 1;
  if (mag_magfield(k) == 3)
    mag_magfield_multiplier = -1;
  end

  conl_ov_mfp(k)   = fulltransp.conl_over_mfp;
  alambda_bb(k)    = mag_magfield_multiplier*efinal.alambda_bb;
  epseff(k)        = mag_magfield_multiplier*efinal.epseff3_2;
  check_onsager(k) = abs(efinal.qflux_e - efinal.qcurr_g)./efinal.qcurr_g;
end

[conl_ov_mfp, sortidx] = sort(conl_ov_mfp);
alambda_bb             = alambda_bb(sortidx);
epseff                 = epseff(sortidx);

figure(1)
loglog(conl_ov_mfp, (epseff), 'o-')
hold on
xlabel('Lc/lc')
ylabel('epseff_3/2', 'Interpreter', 'none')

figure(2)
semilogx(conl_ov_mfp, alambda_bb, 'o-')
hold on
xlabel('Lc/lc')
ylabel('alambda_bb', 'Interpreter', 'none')

figure(3)
semilogx(conl_ov_mfp, check_onsager, '-o')
hold on
xlabel('Lc/lc')
ylabel('abs(qflux_e - qcurr_g)./qcurr_g', 'Interpreter', 'none')
