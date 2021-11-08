clear all;
addpath('~/Programs/neo-2/OctaveScripts/')
prop_list = 5:40;

fp_base = 'fieldpropagator';
fr_base = 'fieldripple';

ext = 'h5';
f_name = ['magnetics','.',ext];
if ~exist('mag_s','var')
  %~ mag_s = hdf52struct(f_name);
  mag_s = h52struct(f_name);
  %~ mag_s = load(f_name);
end

props = fieldnames(mag_s.('fieldpropagator'));
props_num = cellfun(@(x) str2num(strrep(x,'fieldpropagator_','')),props);
%~ props_num = cellfun(@(x) str2num(strrep(x,'_','')),props);
[props_num,s_ind] = sort(props_num);

% for line plots
phi_line = [];
eta_line = [];
phi_bhat_line = [];
bhat_line = [];
phi_bmin_line = [];
bmin_line = [];
phi_bmax_line = [];
bmax_line = [];
phi_ripple_line = [];
ripple_line = [];
phi_prop_line = [];
prop_line = [];

% figure;
count = 0;
ripple_tag_prev = 0;
for prop = prop_list
  if prop<min(props_num), continue, end
  if prop>max(props_num), break; end
  count = count + 1;
  fp_name = ['fieldpropagator_',sprintf('%i',prop)];
  %~ fp_name = ['_',sprintf('%i',prop)];
  phi_l = mag_s.(fp_base).(fp_name).phi_l;
  if count==1, phi_start = phi_l; end
  phi_r = mag_s.(fp_base).(fp_name).phi_r;
  phi = mag_s.(fp_base).(fp_name).x2;
  has_min = mag_s.(fp_base).(fp_name).has_min;

  bhat = mag_s.(fp_base).(fp_name).bhat;
  ripple_tag = mag_s.(fp_base).(fp_name).ch_tag;
  if ripple_tag ~= ripple_tag_prev
    fr_name = ['fieldripple_',sprintf('%i',ripple_tag)];
    %~ fr_name = ['_',sprintf('%i',ripple_tag)];
    eta = mag_s.(fr_base).(fr_name).eta;
    eta_cl = mag_s.(fr_base).(fr_name).eta_cl;
    b_max_l = mag_s.(fr_base).(fr_name).b_max_l;
    b_max_r = mag_s.(fr_base).(fr_name).b_max_r;
    phi_max_l = mag_s.(fr_base).(fr_name).phi_l;
    phi_max_r = mag_s.(fr_base).(fr_name).phi_r;
    b_min = mag_s.(fr_base).(fr_name).b_min;
    phi_min = mag_s.(fr_base).(fr_name).phi_min;
  end
  if phi_l == phi_max_l
    ripple_left_boundary = true;
  else
    ripple_left_boundary = false;
  end
  if phi_r == phi_max_r
    ripple_right_boundary = true;
  else
    ripple_right_boundary = false;
  end

  % eta-levels
  eta_plot = eta;

  %~ numel(eta_plot)
  for k_line = 1:numel(eta_plot)
    phi_ind = find( 1./bhat >= eta_plot(k_line) );
    if ~isempty(phi_ind)
      phi_st = phi(phi_ind(1));
      phi_en = phi(phi_ind(end));
      % old Version with many line objects
      % line([phi_st,phi_en],[eta_plot(k_line),eta_plot(k_line)])
      phi_line = [phi_line;phi_st;phi_en;nan];
      eta_line = [eta_line;eta_plot(k_line);eta_plot(k_line);nan];
    end
  end
  % new version with one line object per

  % b-Field
  if count == 1
    phi_bhat_line = [phi_bhat_line;phi];
    bhat_line = [bhat_line;1./bhat];
  else
    phi_bhat_line = [phi_bhat_line;phi(2:end)];
    bhat_line = [bhat_line;1./bhat(2:end)];
  end

  % line(phi,1./bhat,'Color','r')

  % Minimum
  if has_min
    phi_bmin_line = [phi_bmin_line;phi_min;nan];
    bmin_line = [bmin_line;1./b_min;nan];
    %line(phi_min,1./b_min,'Marker','o','MarkerFaceColor','k',...
    %    'LineStyle','none','MarkerEdgeColor','k','MarkerSize',3)
  end
  % Maxima
  if ripple_left_boundary
    phi_bmax_line = [phi_bmax_line;phi_max_l;nan];
    bmax_line = [bmax_line;1./b_max_l;nan];
    %line([phi_max_l],[1./b_max_l],'Marker','o','MarkerFaceColor','r',...
    %    'LineStyle','none','MarkerEdgeColor','r','MarkerSize',3)
  end
  if ripple_right_boundary
    phi_bmax_line = [phi_bmax_line;phi_max_r;nan];
    bmax_line = [bmax_line;1./b_max_r;nan];
    %line([phi_max_r],[1./b_max_r],'Marker','o','MarkerFaceColor','r',...
    %    'LineStyle','none','MarkerEdgeColor','r','MarkerSize',3)
  end

  % Borders
  if ripple_left_boundary
    phi_ripple_line = [phi_ripple_line;phi_l;phi_l;nan];
    ripple_line = [ripple_line;-1;1;nan];
    %line([phi_l,phi_l],ylim,'Color','r');
  elseif ripple_right_boundary
    phi_ripple_line = [phi_ripple_line;phi_r;phi_r;nan];
    ripple_line = [ripple_line;-1;1;nan];
    %line([phi_r,phi_r],ylim,'Color','r');
  end

  if count==1
    phi_prop_line = [phi_prop_line;phi_l;phi_l;nan];
    prop_line = [prop_line;-1;1;nan];
    %line([phi_l,phi_l],ylim,'Color','b'),
  else
    phi_prop_line = [phi_prop_line;phi_r;phi_r;nan];
    prop_line = [prop_line;-1;1;nan];
    %line([phi_r,phi_r],ylim,'Color','b')
  end
  ripple_tag_prev = ripple_tag;

  phi_end = phi_r;
end

% all plots
figure
line(phi_line,eta_line)
line(phi_bhat_line,bhat_line,'Color','r')
line(phi_bmin_line,bmin_line,'Marker','o','MarkerFaceColor','k',...
    'LineStyle','none','MarkerEdgeColor','k','MarkerSize',3)
line(phi_bmax_line,bmax_line,'Marker','o','MarkerFaceColor','r',...
            'LineStyle','none','MarkerEdgeColor','r','MarkerSize',3)
y_lim = ylim;
ripple_line(ripple_line==-1) = y_lim(1);
ripple_line(ripple_line==+1) = y_lim(2);
prop_line(prop_line==-1) = y_lim(1);
prop_line(prop_line==+1) = y_lim(2);
line(phi_prop_line,prop_line,'Color','b')
line(phi_ripple_line,ripple_line,'Color','r')

xlim([phi_start,phi_end])
xlabel('$\varphi_s$', 'Interpreter', 'latex')
ylabel('$1 / \hat B$', 'Interpreter', 'latex')

%% Export to PGF
[phi_vec, phi_vec_uind] = unique(cell2mat(phi_cell));
theta_vec = cell2mat(theta_cell);
theta_vec = theta_vec(phi_vec_uind);
bhat_vec  = cell2mat(bhat_cell);
bhat_vec  = bhat_vec(phi_vec_uind);
[phi_max_vec, phi_max_uind] = unique(cell2mat(phi_max_cell));
[phi_min_vec, phi_min_uind] = unique(cell2mat(phi_min_cell));
b_max_vec   = cell2mat(b_max_cell);
b_max_vec   = b_max_vec(phi_max_uind);
b_min_vec   = cell2mat(b_min_cell);
b_min_vec   = b_min_vec(phi_min_uind);

pp_bhat = spline(phi_vec, bhat_vec);

plot(phi_vec, bhat_vec)
hold on
plot(phi_min_vec, b_min_vec, 'x')
plot(phi_max_vec, b_max_vec, 'x')

max_ind_oi = 7;
plot(phi_max_vec(max_ind_oi), b_max_vec(max_ind_oi), 'o')
plot(phi_max_vec(max_ind_oi+1), b_max_vec(max_ind_oi+1), 'o')

n_points   = 200*10;
phi_oi_start = phi_max_vec(max_ind_oi)+2*pi*50/mag_s.surface.tag_1.aiota;
phi_oi_end   = phi_oi_start + 2*pi*30/mag_s.surface.tag_1.aiota
phi_b_oi   = linspace(phi_oi_start, phi_oi_end, n_points);
theta_b_oi = theta_vec(1) + (phi_b_oi - phi_vec(1)) * mag_s.surface.tag_1.aiota;

colors = jet(numel(phi_b_oi));
for k = 1:numel(phi_b_oi)
  plot(phi_b_oi(k), ppval(pp_bhat, phi_b_oi(k)), 'o', 'MarkerFaceColor', colors(k,:), 'MarkerSize', 10)
end
hold off

disp('Input-File')
fprintf('%19s %18s %18s %5s\n ', 's', 'theta', 'phi', 'tag')
str = '';
for k = 1:numel(phi_b_oi)
  fprintf('%18.8e %18.8e %18.8e %5s\n ', 0.25d0, theta_b_oi(k), mod(phi_b_oi(k), 2*pi/5), sprintf('p%04d', k))
  str = [str, sprintf('%18.8e %18.8e %18.8e %5s\n ', 0.25d0, theta_b_oi(k), mod(phi_b_oi(k), 2*pi/5), sprintf('p%04d', k))];
end
clipboard('copy', str)

Export to PGF with sparsing
n = 2000;
pp_phi = spline(1:numel(phi_vec), phi_vec);

phi_plot  = interp1(1:numel(phi_vec), phi_vec, linspace(1, numel(phi_vec), n), 'pchip');
bhat_plot = interp1(phi_vec, bhat_vec, phi_plot, 'pchip');

plot_name = 'w7x-sc1-s025-fieldline'
PGFdata = [phi_plot(:) * mag_s.surface.tag_1.aiota, bhat_plot(:)];
save([plot_name, '_pgf.dat'], 'PGFdata', '-ascii')
PGFdata = [phi_b_oi(:).*  mag_s.surface.tag_1.aiota, ppval(pp_bhat, phi_b_oi(:))];
save([plot_name, '_poi_pgf.dat'], 'PGFdata', '-ascii')
