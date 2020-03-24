% \brief Create a perturbed boozer file for given equlibrium input and displacements.
%
%
% Input:
% ------
% file_base: base filename (i.e. without extension) for the input file
%   with the background/unperturbed field.
%   default value: 'eqdsk_35568_2.68800'
% file_displacement: filename (i.e. with extension) for the input file
%   with the displacements caused by the modes.
%   default value: 'displacement_morepoints.txt'
% amplitudes: array with the amplitudes of the single modes. At the
%   moment two values are expected.
%   default value: [1/100, 1/150]
% phases: array with the phases of the single modes. At the moment two
%   values are expected.
%   default value: [0, pi]
% plot_data: logical, if true, plots of the results are made.
%   default value: false;
% file_out_name: name for the output file, without extension.
%   default value: [file_base,'-pert_kink_tear-n',num2str(n_pert),...
%   '-','phase_kt-',num2str(delta_phase),'pi', '.bc']
function create_asdex_perturb_incompress(file_base, file_displacement, amplitudes, phases, plot_data, file_out_name)
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
  % nargin < 6 checked below, default value can not be calculated at
  % this point.

  file_ext = 'bc';

  % Perturbation field
  m_kink = 2;   % poloidal perturbation number for kink
  m_tear = 3;   % poloidal perturbation number for tearing
  n_pert = 2;   % toroidal perturbation number

  % Multiplier for number of modes
  n_multiplier = 4;

  grid_multiplier_theta = 10;
  grid_multiplier_phi = 10;

  tearing_exponent_factor = 1.0/0.06;

  eps = 0.01;

  %% Code

  % coloumns: rho_tor, rho_tor^2, \Delta r_kink, \Delta r_tear
  data_displacement = load(file_displacement);
  rho_tor_dis = data_displacement(:,1);
  s = data_displacement(:, 2);
  delta_rho22 = data_displacement(:,3);
  delta_rho32 = data_displacement(:,4);



  % Define filename
  file_in = [file_base,'.',file_ext];

  delta_phase=(phases(2)-phases(1))/pi;

  if nargin < 6 || isempty(file_out_name)
    file_out = [file_base,'-pert_kink_tear-n',num2str(n_pert),'-','phase_kt-',num2str(delta_phase),'pi'];
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
  phi_vec   = linspace(0, 2*pi, N_phi);
  [T, P] = meshgrid(theta_vec(1:end-1), phi_vec(1:end-1));

  s0 = zeros(nsurf,1);

  for ns = 1:nsurf
%    disp(['Preprocessing flux surface: ', num2str(ns), '/', num2str(nsurf)]);
    head = data_c(k_start+1:k_start+4);
    spec = data_c(k_start+5:k_start+5 + (m0b+1)*(n0b+1) - 1);

    k_start = k_start+5 + (m0b+1)*(n0b+1) - 1;  %for the next flux surface

    % write the original header for flux surface
    for kl = 1:numel(head)
      if (kl == 3)
        headnum = sscanf(head{3},'%f').';
        s0(ns)    = headnum(1);
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
    bmc_f(ns, :) = data(:,9);
    bms_f(ns, :) = data(:,10);
    Rmc_f(ns, :) = data(:,3);
    Rms_f(ns, :) = data(:,4);
    Zmc_f(ns, :) = data(:,5);
    Zms_f(ns, :) = data(:,6);
  end

%cccccccccccccccccc
  xis_kink = 2*rho_tor_dis.*delta_rho22;
  xis_tear = 2*rho_tor_dis.*(delta_rho32 - delta_rho32(1).*exp(-rho_tor_dis*tearing_exponent_factor));
  spl_kink = spline(s, xis_kink);
  spl_tear = spline(s, xis_tear);
  clear xis_kink xis_tear
  xis_kink = ppval(spl_kink, s0);
  xis_tear = ppval(spl_tear, s0);
  der_xis_kink = ppval(fnder(spl_kink), s0);
  der_xis_tear = ppval(fnder(spl_tear), s0);
  max_der_kink = max(abs(der_xis_kink));
  max_der_tear = max(abs(der_xis_tear));

  % Normalize perturbation such that there is only a touching of fieldlines, but no crossing.
  % eps is supposed to be small, and making sure, that there should should be no real touching either.
  xis_kink = xis_kink/(max_der_kink*(1+eps));
  der_xis_kink = der_xis_kink/(max_der_kink*(1+eps));
  xis_tear = xis_tear/(max_der_tear*(1+eps));
  der_xis_tear = der_xis_tear/(max_der_tear*(1+eps));

  % This is the ratio of amplitudes before multiplication with the given amplitudes.
  max_kink_to_max_tear =  max(xis_kink)/max(xis_tear)

  xis_kink = xis_kink*amplitudes(1);
  der_xis_kink = der_xis_kink*amplitudes(1);
  xis_tear = xis_tear*amplitudes(2);
  der_xis_tear = der_xis_tear*amplitudes(2);

  if plot_data
    figure;
    plot(sqrt(s0), xis_kink, sqrt(s0), xis_tear);
  end

  thp = linspace(0, 2*pi, 1000);

  theold = linspace(0, 2*pi, N_theta);
  phi = linspace(0, 2*pi, N_phi);

  for ip=1:1:N_phi

    for ir=1:1:nsurf
      % theta unperturbed on the equidistant grid of theta perturbed:
      thet = thp;
      thet = thet + sin(m_kink*thp - n_pert*phi(ip) + phases(1))*der_xis_kink(ir)/m_kink;
      thet = thet + sin(m_tear*thp - n_pert*phi(ip) + phases(2))*der_xis_tear(ir)/m_tear;

      % theta perturbed on the equidistant grid of theta unperturbed:
      thenew(ir, :, ip) = ppval(spline(thet, thp), theold);

      % s perturbed on the equidistant grid of theta unperturbed:
      snew(ir, :, ip) = s0(ir) + xis_kink(ir).*cos(m_kink*thenew(ir,:,ip) - n_pert*phi(ip) + phases(1)) ...
                               + xis_tear(ir).*cos(m_tear*thenew(ir,:,ip) - n_pert*phi(ip) + phases(2));
    end
    disp(['snew and thenew for ip = ', num2str(ip), '/', num2str(N_phi)]);
  end

  nharm=size(m_orig,1);

  % Fourier amplitudes for perturbed values of s:
  for iharm=1:1:nharm
    bmc_spl(iharm) = spline(s0,bmc_f(:, iharm));
    bms_spl(iharm) = spline(s0,bms_f(:, iharm));
    Rmc_spl(iharm) = spline(s0,Rmc_f(:, iharm));
    Rms_spl(iharm) = spline(s0,Rms_f(:, iharm));
    Zmc_spl(iharm) = spline(s0,Zmc_f(:, iharm));
    Zms_spl(iharm) = spline(s0,Zms_f(:, iharm));
  end

  all_b0_pert = zeros(nsurf, N_phi-1, N_theta-1);
  all_Z_pert = zeros(nsurf, N_phi-1, N_theta-1);
  all_R_pert = zeros(nsurf, N_phi-1, N_theta-1);

  for ip=1:1:N_phi-1
    for it=1:1:N_theta-1
      for iharm=1:1:nharm
        bmc_spert(:, iharm) = ppval(bmc_spl(iharm), snew(:, it, ip));
        bms_spert(:, iharm) = ppval(bms_spl(iharm), snew(:, it, ip));
        Rmc_spert(:, iharm) = ppval(Rmc_spl(iharm), snew(:, it, ip));
        Rms_spert(:, iharm) = ppval(Rms_spl(iharm), snew(:, it, ip));
        Zmc_spert(:, iharm) = ppval(Zmc_spl(iharm), snew(:, it, ip));
        Zms_spert(:, iharm) = ppval(Zms_spl(iharm), snew(:, it, ip));
      end

      for ir=1:1:nsurf
        all_b0_pert(ir, ip, it) = bmc_spert(ir, :)*cos(m_orig*thenew(ir, it, ip)) + bms_spert(ir, :)*sin(m_orig*thenew(ir, it, ip));
        all_R_pert(ir, ip, it) = Rmc_spert(ir, :)*cos(m_orig*thenew(ir, it, ip)) + Rms_spert(ir, :)*sin(m_orig*thenew(ir, it, ip));
        all_Z_pert(ir, ip, it) = Zmc_spert(ir, :)*cos(m_orig*thenew(ir, it, ip)) + Zms_spert(ir, :)*sin(m_orig*thenew(ir, it, ip));
      end
    end
    disp(['perturbed B, R and Z for ip = ', num2str(ip), '/', num2str(N_phi)]);
  end
%
%cccccccccccccccccc
%

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
