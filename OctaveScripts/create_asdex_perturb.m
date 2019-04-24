function create_asdex_perturb(file_base, file_displacement, amplitudes, phases, plot_data)
  %% Parameters
  % Input files
  if nargin < 1 || isempty(file_base)
    file_base = 'eqdsk_35568_2.68800';
  end
  if nargin < 2 || isempty(file_displacement)
    file_displacement = 'displacement_morepoints.txt';
  end
  if nargin < 3 || isempty(amplitudes)
    amplitudes = [1/100, 1/150];
  end
  if nargin < 4 || isempty(phases)
    phases = [0, pi];
  end
  if nargin < 5 || isempty(plot_data)
    plot_data = false;
  end

  file_ext = 'bc';

  % Perturbation field
  m_kink = 2;   % poloidal perturbation number for kink
  m_tear = 3;   % poloidal perturbation number for tearing
  n_pert = 2;   % toroidal perturbation number

  % Multiplier for number of modes
  n_multiplier = 4;

  grid_multiplier = 20; %10; % 100;

  %% Code

  % coloumns: rho_tor, rho_tor^2, \Delta r_kink, \Delta r_tear
  data_displacement = load(file_displacement);
  rho_tor_dis = data_displacement(:,1);
  delta_rho22 = data_displacement(:,3)*amplitudes(1);
  delta_rho32 = data_displacement(:,4)*amplitudes(2);

  if plot_data
    figure;
    plot(rho_tor_dis, rho_tor_dis.^2+2*rho_tor_dis.*delta_rho22, rho_tor_dis, rho_tor_dis.^2+2*rho_tor_dis.*delta_rho32);
    hold on
    plot(rho_tor_dis, (rho_tor_dis+delta_rho22).^2, rho_tor_dis, (rho_tor_dis+delta_rho32).^2);
    plot(rho_tor_dis, (rho_tor_dis-delta_rho22).^2, rho_tor_dis, (rho_tor_dis-delta_rho32).^2);
    hold off
  end

  % Define filename
  file_in = [file_base,'.',file_ext];

  delta_phase=(phases(2)-phases(1))/pi;
  file_out = [file_base,'-pert_kink_tear-n',num2str(n_pert),'-','phase_kt-',num2str(delta_phase),'pi'];
  file_out = [file_out, '.', file_ext];

  % Read Boozer file
  data_c = read_file_into_cell_array(file_in);

  % find m0b and n0b, the value is in the line after the one that contains 'm0b'.
  for k = 1:numel(data_c)
    tline = data_c{k};
    if strfind(tline,'m0b')
      k_def = k + 1;
      break
    end
  end

  % extract m0b and n0b
  tline = data_c{k_def};
  % strtrim is required, because otherwise strsplit would create an
  % empty first entry, due to the whitespace(s) at the beginning of the
  % line.
  dline = str2double(strsplit(strtrim(tline), {' '}));
  m0b = dline(1);
  n0b = dline(2);
  nsurf = dline(3);

  % Number of grid points
  N_theta = m0b*grid_multiplier+2;
  N_phi   = n_pert*grid_multiplier+2;

  % Modify and write

  % change now n0b according to perturbation
  dline(2) = max(n_pert)*n_multiplier/2;
  dline(2) = dline(2) * 2; % because of negative perturbations
  data_c{k_def} = sprintf('   %d    %d   %d    %d %.6E   %.5f   %.5f', dline);

  % open output file and write first lines
  fid = fopen(file_out,'w');
  for kl = 1:k_def
    fprintf(fid,'%s\n',data_c{kl});
  end

  % go through all flux surfaces

  k_start = k_def;

  all_head_data = zeros(nsurf,6);
  all_b0_theta = zeros(nsurf, N_theta-1);
  all_R_theta = zeros(nsurf, N_theta-1);
  all_Z_theta = zeros(nsurf, N_theta-1);
  all_T = zeros(nsurf, N_theta-1); %!!!


  theta_vec = linspace(0, 2*pi, N_theta);
  phi_vec   = linspace(0, 2*pi, N_phi);
  [T, P] = meshgrid(theta_vec(1:end-1), phi_vec(1:end-1));

  s0 = zeros(nsurf,1);
  delta_rho22_A = zeros(nsurf,1);
  delta_rho32_A = zeros(nsurf,1);

  for ns = 1:nsurf
    disp(['Preprocessing flux surface: ', num2str(ns), '/', num2str(nsurf)]);
    head = data_c(k_start+1:k_start+4);
    spec = data_c(k_start+5:k_start+5 + (m0b+1)*(n0b+1) - 1);

    k_start = k_start+5 + (m0b+1)*(n0b+1) - 1;  %for the next flux surface

    % write the original header for flux surface
    for kl = 1:numel(head)
      if (kl == 3)
        headnum = sscanf(head{3},'%f').';
        s0(ns)    = headnum(1);
        delta_rho22_A(ns) =  spline(rho_tor_dis, delta_rho22, sqrt(s0(ns)));
        delta_rho32_A(ns) =  spline(rho_tor_dis, delta_rho32, sqrt(s0(ns)));

        all_head_data (ns,:) = headnum;
      end
    end

    % extract data from spec
    data = zeros(numel(spec),10);
    for kl = 1:numel(spec)
      data(kl,:) = sscanf(spec{kl},'%f').';
    end
    b = data(:,6);
    b00 = b(1);

    % Here are the original data per flux surface
    m_orig  = data(:,1);
    bmc_orig = data(:,9);
    bms_orig = data(:,10);
    Rmc_orig = data(:,3);
    Rms_orig = data(:,4);
    Zmc_orig = data(:,5);
    Zms_orig = data(:,6);

    all_b0_theta(ns, :) = b_fourier(m_orig, bmc_orig, theta_vec(1:end-1))+b_fourier_sin(m_orig, bms_orig, theta_vec(1:end-1));
    all_R_theta(ns, :) = b_fourier(m_orig, Rmc_orig, theta_vec(1:end-1))+b_fourier_sin(m_orig, Rms_orig, theta_vec(1:end-1));
    all_Z_theta(ns, :) = b_fourier(m_orig, Zmc_orig, theta_vec(1:end-1))+b_fourier_sin(m_orig, Zms_orig, theta_vec(1:end-1));
    all_T(ns,:) = theta_vec(1:end-1); %!!!

  end

  if plot_data
    figure
    plot(s0,delta_rho32_A)
  end

  all_b0_unpert = zeros(nsurf, N_phi-1, N_theta-1);
  all_R_unpert = zeros(nsurf, N_phi-1, N_theta-1);
  all_Z_unpert = zeros(nsurf, N_phi-1, N_theta-1);

  all_b0_pert = zeros(nsurf, N_phi-1, N_theta-1);
  all_Z_pert = zeros(nsurf, N_phi-1, N_theta-1);
  all_R_pert = zeros(nsurf, N_phi-1, N_theta-1);

  all_b0_pert2 = zeros(nsurf, N_phi-1, N_theta-1);

  all_T_tri = zeros(nsurf, N_phi-1, N_theta-1);

  for i_phi = 1: N_phi-1
   all_b0_unpert(:,i_phi,:) = all_b0_theta;
   all_R_unpert(:,i_phi,:) = all_R_theta;
   all_Z_unpert(:,i_phi,:) = all_Z_theta;
   all_T_tri(:,i_phi,:) = all_T;
  end

  cos_var_kink = cos(m_kink*T - n_pert*P + phases(1));
  cos_var_tear = cos(m_tear*T - n_pert*P + phases(2));


  for i_theta = 1: N_theta-1
    [sqrt_so_mx, ~ ] = meshgrid(sqrt(s0),cos_var_kink(:,i_theta));
    [delta_rho22_A_mx, local_cos_var_kink] = meshgrid(delta_rho22_A,cos_var_kink(:,i_theta));
    [delta_rho32_A_mx, local_cos_var_tear] = meshgrid(delta_rho32_A,cos_var_tear(:,i_theta));
    boozer_s_pert = (sqrt_so_mx + delta_rho22_A_mx .* local_cos_var_kink + delta_rho32_A_mx .* local_cos_var_tear).^2;

    % Explicitly set the behaviour, when extrapolation is necessary
    % (i.e. at the borders in radial direction), otherwise matlab and
    % octave might do different things (experienced in 2018b vs. 4.0.3).
    local_b0_pert = interp1(s0, all_b0_unpert(:,1,i_theta), boozer_s_pert, 'spline', 'extrap');
    local_R_pert  = interp1(s0, all_R_unpert(:,1,i_theta),  boozer_s_pert, 'spline', 'extrap');
    local_Z_pert  = interp1(s0, all_Z_unpert(:,1,i_theta),  boozer_s_pert, 'spline', 'extrap');
    all_b0_pert(:,:,i_theta) = local_b0_pert' ;
    all_R_pert(:,:,i_theta)  = local_R_pert' ;
    all_Z_pert(:,:,i_theta)  = local_Z_pert' ;
  end

  T_tri = all_T_tri(1, :, :);
  if sum(T(:)==T_tri(:)) ~= numel(T)
    error('Error: b0_pert has corrupt dimension indices!')
  end

  k_start = k_def;

  for ns = 1:nsurf
    disp(['Processing flux surface: ', num2str(ns), '/', num2str(nsurf)]);
    head = data_c(k_start+1:k_start+4);
    spec = data_c(k_start+5:k_start+5 + (m0b+1)*(n0b+1) - 1);
    k_start = k_start+5 + (m0b+1)*(n0b+1) - 1;  %for the next flux surface

    % write the original header for flux surface
    for kl = 1:numel(head)
      if (kl == 3)
        headnum = all_head_data (ns,:);
        fprintf(fid,' %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n',headnum);
      else
        fprintf(fid,'%s\n',head{kl});
      end
    end

    % extract data from spec
    data = zeros(numel(spec),10);
    for kl = 1:numel(spec)
      data(kl,:) = sscanf(spec{kl},'%f').';
    end
    b = data(:,6);
    b00 = b(1);

    % Here are the original data per flux surface
    m_orig  = data(:,1);
    bm_orig = data(:,6);

    mm_max = m0b;
    nn_max = n_pert*n_multiplier;

    m_vec  = 0:mm_max;%-mm_max:1:mm_max;
    n_vec  = -nn_max:1:nn_max;

    bmnc    = zeros(numel(m_vec), numel(n_vec));
    bmns    = zeros(numel(m_vec), numel(n_vec));
    Rmnc    = zeros(numel(m_vec), numel(n_vec));
    Rmns    = zeros(numel(m_vec), numel(n_vec));
    Zmnc    = zeros(numel(m_vec), numel(n_vec));
    Zmns    = zeros(numel(m_vec), numel(n_vec));

    b0_pert = all_b0_pert(ns, :, :);
    R_pert  = all_R_pert(ns, :, :);
    Z_pert  = all_Z_pert(ns, :, :);
    % Compute bmn
    for mm = 1:numel(m_vec)
      for nn = 1:numel(n_vec)
        m = m_vec(mm);
        n = n_vec(nn);

        if (mod(n, n_pert) == 0)
          bmnc(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(cos(m*T(:) - n*P(:)) .* b0_pert(:));
          bmns(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(sin(m*T(:) - n*P(:)) .* b0_pert(:));
          Rmnc(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(cos(m*T(:) - n*P(:)) .* R_pert(:));
          Rmns(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(sin(m*T(:) - n*P(:)) .* R_pert(:));
          Zmnc(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(cos(m*T(:) - n*P(:)) .* Z_pert(:));
          Zmns(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(sin(m*T(:) - n*P(:)) .* Z_pert(:));
        else
          bmnc(mm, nn) = 0;
          bmns(mm, nn) = 0;
          Rmnc(mm, nn) = 0;
          Rmns(mm, nn) = 0;
          Zmnc(mm, nn) = 0;
          Zmns(mm, nn) = 0;
        end
      end
    end

    % Write
    data_orig = data;
    for n_p = -nn_max:nn_max
      data = data_orig;
      data(:,2) = n_p;

      data(:,3) = Rmnc(:, n_p+nn_max+1);
      data(2:end,3) = 2*data(2:end,3);
      data(:,4) = 2*Rmns(:, n_p+nn_max+1);

      data(:,5) = Zmnc(:, n_p+nn_max+1);
      data(2:end,5) = 2*data(2:end,5);
      data(:,6) = 2*Zmns(:, n_p+nn_max+1);


      data(:,9) = bmnc(:, n_p+nn_max+1);
      data(2:end,9) = 2*data(2:end,9);
      data(:,10) = 2*bmns(:, n_p+nn_max+1);


      if (n_p ~= 0)
        data(:,7:8) = 0;
      end

      for kl = 1:size(data,1)
        fprintf(fid,'   %2d   %2d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n', data(kl,:));
      end
    end
  end
  fclose(fid);
end

function [data] = read_file_into_cell_array(filename)
  fid = fopen(filename);

  k = 1;
  tline = fgetl(fid);
  data{k} = tline;
  while ischar(tline)
    tline = fgetl(fid);
    if ischar(tline)
      k = k + 1;
      data{k} = tline;
    end
  end

  fclose(fid);
end
