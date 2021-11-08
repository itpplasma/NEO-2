%matlab -nodesktop -nosplash

%%Preloader
%close all

addpath('/proj/plasma/Neo2/Interface/Matlab/')

fp_base = 'fieldpropagator';
fr_base = 'fieldripple';

ext = 'h5';
f_name = ['','magnetics','.',ext];
if ~exist('mag_s','var')
  mag_s = hdf52struct(f_name);
end

%%
props = fieldnames(mag_s.('fieldpropagator'));
props_num = cellfun(@(x) str2num(strrep(x,'tag_','')),props);
[props_num,s_ind] = sort(props_num);
prop_list = reshape(props_num,1,[]);

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

count = 0;
ripple_tag_prev = 0;
for prop = prop_list(1:end)
  if prop<min(props_num), continue, end
  if prop>max(props_num), break; end
  count = count + 1;
  fp_name = ['tag_',sprintf('%i',prop)];
  phi_l = mag_s.(fp_base).(fp_name).phi_l;
  if count==1, phi_start = phi_l; end
  phi_r = mag_s.(fp_base).(fp_name).phi_r;
  phi = mag_s.(fp_base).(fp_name).x2;
  has_min = mag_s.(fp_base).(fp_name).has_min;

  bhat = mag_s.(fp_base).(fp_name).bhat;
  ripple_tag = mag_s.(fp_base).(fp_name).ch_tag;
  if ripple_tag ~= ripple_tag_prev
    fr_name = ['tag_',sprintf('%i',ripple_tag)];
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

  % b-Field
  if count == 1
    phi_bhat_line = [phi_bhat_line;phi];
    bhat_line = [bhat_line;1./bhat];
  else
    phi_bhat_line = [phi_bhat_line;phi(2:end)];
    bhat_line = [bhat_line;1./bhat(2:end)];
  end

  % Minumum
  if has_min
    phi_bmin_line = [phi_bmin_line;phi_min];
    bmin_line = [bmin_line;b_min];
    %line(phi_min,1./b_min,'Marker','o','MarkerFaceColor','k',...
    %    'LineStyle','none','MarkerEdgeColor','k','MarkerSize',3)
  end
  % Maxima
  if ripple_left_boundary
    phi_bmax_line = [phi_bmax_line;phi_max_l];
    bmax_line = [bmax_line;b_max_l];
    %line([phi_max_l],[1./b_max_l],'Marker','o','MarkerFaceColor','r',...
    %    'LineStyle','none','MarkerEdgeColor','r','MarkerSize',3)
  end
  if ripple_right_boundary
    phi_bmax_line = [phi_bmax_line;phi_max_r];
    bmax_line = [bmax_line;b_max_r];
    %line([phi_max_r],[1./b_max_r],'Marker','o','MarkerFaceColor','r',...
    %    'LineStyle','none','MarkerEdgeColor','r','MarkerSize',3)
  end

end

numel(phi_bmax_line)
numel(bmax_line)

phi_bmax_line_orig = phi_bmax_line;
bmax_line_orig = bmax_line;

%%
figure
line(phi_bhat_line,bhat_line,'Color','r')
export_pgf(gca, 'bhat_pgf.dat', 1, true)

%%
bmax_line = bmax_line_orig;
num_max1   = numel(bmax_line);

bmax_line_save = bmax_line;

phi_bmax_line_save = phi_bmax_line;

nc_right = cell(1, numel(bmax_line));
kc_right = cell(1, numel(bmax_line));
mc_right = cell(1, numel(bmax_line));
pc_right = cell(1, numel(bmax_line));

nc_left  = cell(1, numel(bmax_line));
kc_left  = cell(1, numel(bmax_line));
mc_left  = cell(1, numel(bmax_line));
pc_left  = cell(1, numel(bmax_line));

[bmax, bmax_ind] = max(bmax_line);

left  = true;
right = true;
bmax_vec = 1:numel(bmax_line_save);
for bmax_ind = bmax_vec

  bmax = bmax_line(bmax_ind);
  num_max1   = numel(bmax_line);
  bmax_line = bmax_line_save;
  bmax_before = bmax_line(bmax_ind+1:end);
  bmax_after  = bmax_line(1:bmax_ind-1);
  bmax_line = [bmax_before; bmax_line; bmax_after];
  bmax_ind_for = numel(bmax_before) + bmax_ind;

  phi_bmax_line = phi_bmax_line_save;
  phi_bmax_before = phi_bmax_line(bmax_ind+1:end);
  phi_bmax_after  = phi_bmax_line(1:bmax_ind-1);
  phi_bmax_line = [phi_bmax_before; phi_bmax_line; phi_bmax_after];

  if (right)
    n = [];
    k = [];
    m = [];
    p = [];
    c1 = 0;
    c2 = 0;
    bmax_block = 0;
    for ni = [bmax_ind_for+1:numel(bmax_line)]
      c1 = c1 + 1;
      if (bmax_line(ni) > bmax_block)
        bmax_block = bmax_line(ni);
        c2 = c2 + 1;
        n(c2) = c1;
        k(c2) = c2;
        m(c2) = bmax_line(ni);
        p(c2) = phi_bmax_line(ni);
      end
    end
    nc_right{bmax_ind} = n;
    kc_right{bmax_ind} = k;
    mc_right{bmax_ind} = m;
    pc_right{bmax_ind} = p;
  end

  if (left)
    n = [];
    k = [];
    m = [];
    p = [];
    c1 = 0;
    c2 = 0;
    bmax_block = 0;
    for ni = [bmax_ind_for-1:-1:1]
      c1 = c1 + 1;
      if (bmax_line(ni) > bmax_block)
        bmax_block = bmax_line(ni);
        c2 = c2 + 1;
        n(c2) = c1;
        k(c2) = c2;
        m(c2) = bmax_line(ni);
        p(c2) = phi_bmax_line(ni);
      end
    end
    nc_left{bmax_ind} = n;
    kc_left{bmax_ind} = k;
    mc_left{bmax_ind} = m;
    pc_left{bmax_ind} = p;
  end
end

phi_bmax_line = phi_bmax_line_save;

%% Observation points
%point = 3800;

% Please note that changing the observation point requires to run again
% the interface!! Otherwise old data will be plotted!
point_phimfl = 4;

point = find(phi_bhat_line - point_phimfl > 0, 1);
phi_l_point  = phi_bhat_line(point);
bhat_point = bhat_line(point);

hold on
plot(phi_l_point, bhat_point, '.k', 'MarkerSize', 30)
hold off
export_pgf(gca, 'poi_pgf.dat', 1)

iota    = mag_s.surface.tag_1.aiota;
theta_0 = mag_s.fieldline.tag_1.xstart(3);
phi_0   = mag_s.fieldline.tag_1.xstart(2);

theta_point = theta_0 + iota * (phi_l_point - phi_0);
phi_point   = phi_l_point;

theta_point = mod(theta_point, 2*pi);
phi_point   = phi_point;%mod(phi_point, 2*pi)

phi_idx_cl = find(diff(sign(phi_bmax_line - phi_l_point)) > 0);
phi_idx_cl_l = phi_idx_cl - 1;
phi_idx_cl_r = phi_idx_cl + 1;

%hold on
%plot(phi_bmax_line(phi_idx_cl_l), 1./bmax_line_orig(phi_idx_cl_l), 'o')
%plot(phi_bmax_line(phi_idx_cl_r), 1./bmax_line_orig(phi_idx_cl_r), 'o')
hold on

pc_vec_x = [pc_right{phi_idx_cl_l}, pc_left{phi_idx_cl_l}, pc_right{phi_idx_cl_r}, pc_left{phi_idx_cl_r}];
pc_vec_y = [1./mc_right{phi_idx_cl_l}, 1./mc_left{phi_idx_cl_l}, 1./mc_right{phi_idx_cl_r}, 1./mc_left{phi_idx_cl_r}];

[pc_vec_x, un_idx] = unique(pc_vec_x);
pc_vec_y = pc_vec_y(un_idx);

%plot(pc_right{phi_idx_cl_l}, 1./mc_right{phi_idx_cl_l}, 'x', 'MarkerSize', 10)
%plot(pc_left{phi_idx_cl_l}, 1./mc_left{phi_idx_cl_l}, '+', 'MarkerSize', 10)
%plot(pc_right{phi_idx_cl_r}, 1./mc_right{phi_idx_cl_r}, '^', 'MarkerSize', 10)
%plot(pc_left{phi_idx_cl_r}, 1./mc_left{phi_idx_cl_r}, 'v', 'MarkerSize', 10)
plot(pc_vec_x, pc_vec_y, 'x', 'MarkerSize', 10)
hold off

export_pgf(gca, 'relmax_pgf.dat', 1);


%%
disp('Input-File')
fprintf('%19s %18s %18s %7s\n ', 's', 'theta', 'phi', 'tag')
str = '';
fprintf('%18.8e %18.8e %18.8e %7s\n ', 0.25d0, theta_point, mod(phi_point, 2*pi), sprintf('p%d', 1))
clipboard('copy', str)


%%
plot_path = [pwd, '/PLOTS/'];
all_rel_max = [1./mc_right{phi_idx_cl_l}, 1./mc_right{phi_idx_cl_r}, 1./mc_left{phi_idx_cl_l}, 1./mc_left{phi_idx_cl_r}];

figure
try
  func_plot_dentf({plot_path}, 1, all_rel_max);
catch
  error('Failed to plot distribution function. Did you run the interface?')
end
hold on


%%
g = gca;
export_pgf(g.Children(1), 'dentf_relmax_pgf.dat')
export_pgf(g.Children(2), 'dentf_tp_pgf.dat')
export_pgf(g.Children(3), 'dentf_dist_pgf.dat', [], true, [0.9, 1])


lh = [];
ls = cell(0);
l  = 0;
for k = 1:numel(g.Children)
  if (~isempty(g.Children(k).DisplayName))
    l = l + 1;
    lh(l) = g.Children(k);
    ls{l} = g.Children(k).DisplayName;
  end
end
legend(lh, ls)

%%
% iplot_vec = 1;
% left = false;
% right = true;
% for iplot = iplot_vec
%
%   %iplot = 3668;%2137; %bmax_ind;%numel(bmax_line); %bmax_ind;
%   figure
%   axes('XScale', 'log')
%   if (right)
%     stairs(nc_right{iplot}, kc_right{iplot}, 'r-x')
%     set(gca, 'XScale', 'log')
%     hold on
%     semilogx(factorial(kc_right{iplot}+1), kc_right{iplot}, 'r-')
%     semilogx(nc_right{iplot}, log(nc_right{iplot}), 'k--')
%   end
%
%   if (left)
%     stairs(nc_left{iplot}, kc_left{iplot}, 'b-o')
%     set(gca, 'XScale', 'log')
%     hold on
%     semilogx(factorial(kc_left{iplot}+1), kc_left{iplot}, 'b:')
%     semilogx(factorial(kc_right{iplot}+1), kc_right{iplot}, 'r-')
%   end
%   xlim([0, max([nc_left{iplot}, nc_right{iplot}])])
%   ylim([0, max([kc_left{iplot}, kc_right{iplot}])])
%   hold off
%   xlabel('n')
%   ylabel('k')
%   legend('Field line', 'Analytical', 'Location', 'Best')
%
% %   figure
% %   plot(bmax_line_orig,'ob-')
% %   hold on
% %   plot(iplot, bmax_line_orig(iplot), 'r*')
% %   hold off
% %   xlim([iplot - 100, iplot + 100])
%
%   if (left)
%     PGFData = [nc_left{iplot}(:), kc_left{iplot}(:)];
%     %save(['rel_maxima_steps_left_', num2str(iplot), '_pgf.dat'], 'PGFData', '-ascii')
%   end
%   if (right)
%     PGFData = [nc_right{iplot}(:), kc_right{iplot}(:)];
%     %save(['rel_maxima_steps_right_', num2str(iplot), '_pgf.dat'], 'PGFData', '-ascii')
%   end
% end
