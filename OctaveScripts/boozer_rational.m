% Script for creating a boozer file with toroidal
% perturbation.

%% Parameters
% Input file
file_base = 'tok-synch2';
file_ext = 'bc';

% Perturbation field
n_pert = 10;  % toroidal perturbation
A_pert = 0.2; % perturbation amplitude

% Multiplier for number of modes
n_multiplier = 4;

% Flat iota profile
lsw_flatiota = true;
flatiota = 4/9;

% Number of grid points
N_theta = n_pert*101;
N_phi   = n_pert*100;

%% Code

% Define filename
file_in = [file_base,'.',file_ext];
s_pert = sprintf('%.0e',A_pert);
file_out = [file_base,'-',num2str(n_pert),'-',strrep(s_pert,'e-0','m'),'-rat'];

if (lsw_flatiota)
  file_out = [file_out, '-iota-', strrep(strtrim(rats(flatiota)), '/', 'ov'), '.', file_ext];
else
  file_out = [file_out, '.', file_ext];
end

%file_out = 'test.bc';

% Read Boozer file
fid = fopen(file_in);

k = 1;
tline = fgetl(fid);
data_c{k} = tline;
while ischar(tline)
  tline = fgetl(fid);
  if ischar(tline)
    k = k + 1;
    data_c{k} = tline;
  end
end
fclose(fid);

% find m0b and n0b
for k = 1:numel(data_c)
  tline = data_c{k};
  if strfind(tline,'m0b')
    k_def = k + 1;
    break
  end
end

% extract m0b and n0b
tline = data_c{k_def};
dline = str2num(tline);
m0b = dline(1);
n0b = dline(2);
nsurf = dline(3);


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
for ns = 1:nsurf
  disp(['Processing flux surface: ', num2str(ns), '/', num2str(nsurf)]);
  head = data_c(k_start+1:k_start+4);
  spec = data_c(k_start+5:k_start+5 + (m0b+1)*(n0b+1) - 1);
  k_start = k_start+5 + (m0b+1)*(n0b+1) - 1;

  % write the original header for flux surface
  for kl = 1:numel(head)
    if (kl == 3) && (lsw_flatiota)
      headnum = sscanf(head{3},'%f').';
      headnum(2) = flatiota;
      fprintf(fid,' %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E\n',headnum);
    else
      fprintf(fid,'%s\n',head{kl});
    end
  end

  % extract data from spec
  data = zeros(numel(spec),6);
  for kl = 1:numel(spec)
    data(kl,:) = sscanf(spec{kl},'%f').';
  end
  b = data(:,6);
  b00 = b(1);

  % Here are the original data per flux surface
  m_orig  = data(:,1);
  bm_orig = data(:,6);

  theta_vec = linspace(0, 2*pi, N_theta);
  phi_vec   = linspace(0, 2*pi, N_phi);
  [T, P] = meshgrid(theta_vec(1:end-1), phi_vec(1:end-1));

  headnum = sscanf(head{3},'%f').';
  s    = headnum(1);
  iota = headnum(2);

  b0_theta = b_fourier(m_orig, bm_orig, T(:));
  b0_pert  = b0_theta ./ sqrt(1 + A_pert * (b0_theta./b00).^2 .* cos(n_pert.* P(:)));

  % This is only for plotting a specific flux surface at the end of the
  % script
  if (abs(s - 2.5781E-01) < 0.001)
    iota_oi = iota;
    m_orig_oi  = m_orig;
    bm_orig_oi = bm_orig;
    b0_theta_oi = b0_theta;
    b0_pert_oi  = b0_pert;
    b00_oi      = b00;
  end

  mm_max = m0b;
  nn_max = n_pert*n_multiplier;

  m_vec  = 0:mm_max;%-mm_max:1:mm_max;
  n_vec  = -nn_max:1:nn_max;

  bmn    = zeros(numel(m_vec), numel(n_vec));

  % Compute bmn
  for mm = 1:numel(m_vec)
    for nn = 1:numel(n_vec)
      m = m_vec(mm);
      n = n_vec(nn);

      if (mod(n, n_pert) == 0)
        bmn(mm, nn) = 1/(N_theta-1) * 1/(N_phi-1) * sum(cos(m*T(:) - n*P(:)) .* b0_pert(:));
      else
        bmn(mm, nn) = 0;
      end
    end
  end

  % Write
  data_orig = data;
  for n_p = -nn_max:nn_max
    data = data_orig;
    data(:,2) = n_p;

    data(:,6) = bmn(:, n_p+nn_max+1);
    data(2:end,6) = 2*data(2:end,6);

    if (n_p ~= 0)
      data(:,3:5) = 0;
    end

    for kl = 1:size(data,1)
      fprintf(fid,'   %2d   %2d %15.8E %15.8E %15.8E %15.8E\n', data(kl,:));
    end
  end
end

fclose(fid);

%% Plotting (only for debugging)
theta = linspace(0, 40*pi, 1000);
phi   = theta / iota_oi;
b0 = b_fourier(m_orig_oi, bm_orig_oi, theta);
b  = b0 ./ sqrt(1 + A_pert * (b0 / b00_oi).^2 .* cos(n_pert * phi));

surf(P, T, reshape(b0_pert_oi, size(P)), 'EdgeColor', 'none')
xlabel('phi');
ylabel('theta');
view(2);

%plot(phi, b, '-');
%hold on

%figure
%plot(theta, 1./b * b00);
