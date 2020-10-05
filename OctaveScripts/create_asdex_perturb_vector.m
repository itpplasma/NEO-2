% \brief Create a perturbed (tearing) boozer file for given equlibrium input.
%
%
% Input:
% ------
% file_base: base filename (i.e. without extension) for the input file
%   with the background/unperturbed field.
%   default value: 'eqdsk_35568_2.68800'
% plot_data: logical, if true, plots of the results are made.
%   default value: false;
% file_out_name: name for the output file, without extension.
%   default value: [file_base,'-pert_kink_tear-n',num2str(n_pert),...
%   '-','phase_kt-',num2str(delta_phase),'pi', '.bc']
function create_asdex_perturb_vector(file_base, plot_data, file_out_name)

  %% Parameters
  % Input files
  if nargin < 1 || isempty(file_base)
    file_base = 'eqdsk_35568_2.68800';
  end
  if nargin < 2 || isempty(plot_data)
    plot_data = false;
  end
  % nargin < 3 checked below, default value can not be calculated at
  % this point.

  file_ext = 'bc';

  % Perturbation field
  m_kink = 2;   % poloidal perturbation number for kink
  m_tear = 3;   % poloidal perturbation number for tearing
  n_pert = 2;   % toroidal perturbation number

  % Multiplier for number of modes
  n_multiplier = 8 %4;

  grid_multiplier_phi = 50; %10;
  grid_multiplier_theta = 10;

  %% Code

  % Define filename
  file_in = [file_base,'.',file_ext];

  if nargin < 3 || isempty(file_out_name)
    file_out = [file_base,'-pert_kink_tear-n',num2str(n_pert)];
    file_out = [file_out, '.', file_ext];
  else
    file_out = file_out_name;
  end

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
  N_theta = m0b*grid_multiplier_theta+2;
  N_phi   = n_pert*grid_multiplier_phi+2;

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


  theta_vec = linspace(0, 2*pi, N_theta);
  phi_vec = linspace(0, 2*pi, N_phi);
  [T, P] = meshgrid(theta_vec(1:end-1), phi_vec(1:end-1));

  s0 = zeros(nsurf,1);
  iota0 = zeros(nsurf,1);

  for ns = 1:nsurf
    % disp(['Preprocessing flux surface: ', num2str(ns), '/', num2str(nsurf)]);
    head = data_c(k_start+1:k_start+4);
    spec = data_c(k_start+5:k_start+5 + (m0b+1)*(n0b+1) - 1);

    k_start = k_start+5 + (m0b+1)*(n0b+1) - 1;  %for the next flux surface

    % write the original header for flux surface
    for kl = 1:numel(head)
      if (kl == 3)
        headnum = sscanf(head{3},'%f').';
        s0(ns) = headnum(1);
        iota0(ns) = headnum(2);
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

    bmc_f(ns,:) = bmc_orig;
    bms_f(ns,:) = bms_orig;
    Rmc_f(ns,:) = Rmc_orig;
    Rms_f(ns,:) = Rms_orig;
    Zmc_f(ns,:) = Zmc_orig;
    Zms_f(ns,:) = Zms_orig;

  end

  % Create perturbed s and theta for tearing mode starting from vector potential:
  vpt

  snew = zeros(nsurf,N_theta-1,N_phi-1);
  thenew = zeros(nsurf,N_theta-1,N_phi-1);
  for ir=1:1:nsurf
    s_and_theta_on_angular_grid
    snew(ir,:,:) = s_surf';
    thenew(ir,:,:) = theta_surf';
    % modify iota profile:
    all_head_data(ir,2) = iota_resc(ir);
  end
  % End create perturbed s and theta for tearing mode starting from vector potential:

  nharm = size(m_orig,1);

  % Fourier amplitudes for perturbed values of s:
  for iharm=1:1:nharm
    bmc_spl(iharm) = spline(s0,bmc_f(:,iharm));
    bms_spl(iharm) = spline(s0,bms_f(:,iharm));
    Rmc_spl(iharm) = spline(s0,Rmc_f(:,iharm));
    Rms_spl(iharm) = spline(s0,Rms_f(:,iharm));
    Zmc_spl(iharm) = spline(s0,Zmc_f(:,iharm));
    Zms_spl(iharm) = spline(s0,Zms_f(:,iharm));
  end

  all_b0_pert = zeros(nsurf, N_phi-1, N_theta-1);
  all_Z_pert = zeros(nsurf, N_phi-1, N_theta-1);
  all_R_pert = zeros(nsurf, N_phi-1, N_theta-1);

  for ip=1:1:N_phi-1
    for it=1:1:N_theta-1
      for iharm=1:1:nharm
        bmc_spert(:,iharm) = ppval(bmc_spl(iharm),snew(:,it,ip));
        bms_spert(:,iharm) = ppval(bms_spl(iharm),snew(:,it,ip));
        Rmc_spert(:,iharm) = ppval(Rmc_spl(iharm),snew(:,it,ip));
        Rms_spert(:,iharm) = ppval(Rms_spl(iharm),snew(:,it,ip));
        Zmc_spert(:,iharm) = ppval(Zmc_spl(iharm),snew(:,it,ip));
        Zms_spert(:,iharm) = ppval(Zms_spl(iharm),snew(:,it,ip));
      end

      for ir=1:1:nsurf
        all_b0_pert(ir,ip,it) = reshape(bmc_spert(ir,:),1,nharm)*cos(m_orig*thenew(ir,it,ip))+reshape(bms_spert(ir,:),1,nharm)*sin(m_orig*thenew(ir,it,ip));
        all_R_pert(ir,ip,it) = reshape(Rmc_spert(ir,:),1,nharm)*cos(m_orig*thenew(ir,it,ip))+reshape(Rms_spert(ir,:),1,nharm)*sin(m_orig*thenew(ir,it,ip));
        all_Z_pert(ir,ip,it) = reshape(Zmc_spert(ir,:),1,nharm)*cos(m_orig*thenew(ir,it,ip))+reshape(Zms_spert(ir,:),1,nharm)*sin(m_orig*thenew(ir,it,ip));
      end
    end
    disp(['perturbed B, R and Z for ip = ', num2str(ip), '/', num2str(N_phi)]);
  end

  k_start = k_def;

  for ns = 1:nsurf
    if(mod(ns,100) == 0)
      disp(['Processing flux surface: ', num2str(ns), '/', num2str(nsurf)]);
    end
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

    bmnc = zeros(numel(m_vec), numel(n_vec));
    bmns = zeros(numel(m_vec), numel(n_vec));
    Rmnc = zeros(numel(m_vec), numel(n_vec));
    Rmns = zeros(numel(m_vec), numel(n_vec));
    Zmnc = zeros(numel(m_vec), numel(n_vec));
    Zmns = zeros(numel(m_vec), numel(n_vec));

    b0_pert = reshape(all_b0_pert(ns, :, :),N_phi-1,N_theta-1);
    R_pert  = all_R_pert(ns, :, :);
    Z_pert  = all_Z_pert(ns, :, :);
    % Compute bmn
    for mm = 1:numel(m_vec)
      for nn = 1:numel(n_vec)
        m = m_vec(mm);
        n = n_vec(nn);

        if (mod(n, n_pert) == 0)
          bmnc(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(cos(m*T(:) + n*P(:)) .* b0_pert(:));
          bmns(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(sin(m*T(:) + n*P(:)) .* b0_pert(:));
          Rmnc(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(cos(m*T(:) + n*P(:)) .* R_pert(:));
          Rmns(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(sin(m*T(:) + n*P(:)) .* R_pert(:));
          Zmnc(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(cos(m*T(:) + n*P(:)) .* Z_pert(:));
          Zmns(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(sin(m*T(:) + n*P(:)) .* Z_pert(:));
          if(m == 0)
            if(n < 0)
              bmnc(mm, nn) = 0;
              bmns(mm, nn) = 0;
              Rmnc(mm, nn) = 0;
              Rmns(mm, nn) = 0;
              Zmnc(mm, nn) = 0;
              Zmns(mm, nn) = 0;
            end
          end
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
