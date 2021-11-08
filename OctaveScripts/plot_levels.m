addpath('/proj/plasma/Neo2/Interface/Matlab/')

%clear all;

% Two periods
prop_list = 2:10000;

% One periodeta
plot_suf = '';
export_as_pgf = false;

%devicename = 'w7x-sc1-s025';
%prop_list = 5:40;

devicename = 'test';
%prop_list = 1:20;

fp_base = 'fieldpropagator';
fr_base = 'fieldripple';

ext = 'h5';
f_name = ['magnetics','.',ext];
%if ~exist('mag_s','var')
    mag_s = hdf52struct(f_name);
%end

props = fieldnames(mag_s.('fieldpropagator'));
props_num = cellfun(@(x) str2num(strrep(x,'tag_','')),props);
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
eta_num = [];
min_diff_eta = [];
% figure;
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
    %numel(eta)
    eta_num(end+1) = numel(eta);
    diff_eta = diff(eta);
    min_diff_eta(end+1) = min(diff_eta);

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
  %

  % b-Field
  if count == 1
    phi_bhat_line = [phi_bhat_line;phi];
    bhat_line = [bhat_line;1./bhat];
  else
    phi_bhat_line = [phi_bhat_line;phi(2:end)];
    bhat_line = [bhat_line;1./bhat(2:end)];
  end



  % line(phi,1./bhat,'Color','r')

  % Minumum
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

%% All plots
figure
line(phi_line,eta_line, 'AlignVertexCenters', 'on')

if (export_as_pgf)
  plot_name = [devicename, '-phi_line'];
  PGFdata = [phi_line(:) - phi_line(1), eta_line(:)];
  save([plot_name, plot_suf, '_pgf.dat'], 'PGFdata', '-ascii')
end

line(phi_bhat_line,bhat_line,'Color','r')
if (export_as_pgf)
  plot_name = [devicename, '-phi_bhat_line'];
  PGFdata = [phi_bhat_line(:) - phi_line(1), bhat_line(:)];
  save([plot_name, plot_suf,'_pgf.dat'], 'PGFdata', '-ascii')
end

line(phi_bmin_line,bmin_line,'Marker','o','MarkerFaceColor','k',...
    'LineStyle','none','MarkerEdgeColor','k','MarkerSize',3)
if (export_as_pgf)
  plot_name = [devicename, '-phi_bmin_line'];
  PGFdata = [phi_bmin_line(:) - phi_line(1), bmin_line(:)];
  save([plot_name, plot_suf,'_pgf.dat'], 'PGFdata', '-ascii')
end

line(phi_bmax_line,bmax_line,'Marker','o','MarkerFaceColor','r',...
    'LineStyle','none','MarkerEdgeColor','r','MarkerSize',3)
if (export_as_pgf)
  plot_name = [devicename, '-phi_bmax_line'];
  PGFdata = [phi_bmax_line(:) - phi_line(1), bmax_line(:)];
  save([plot_name, plot_suf,'_pgf.dat'], 'PGFdata', '-ascii')
end

y_lim = ylim;
ripple_line(ripple_line==-1) = y_lim(1);
ripple_line(ripple_line==+1) = y_lim(2);
prop_line(prop_line==-1) = y_lim(1);
prop_line(prop_line==+1) = y_lim(2);
line(phi_prop_line,prop_line,'Color','b')
plot_name = [devicename, '-phi_prop_line'];
if (export_as_pgf)
  PGFdata = [phi_prop_line(:) - phi_line(1), prop_line(:)];
  save([plot_name, plot_suf,'_pgf.dat'], 'PGFdata', '-ascii')
end

line(phi_ripple_line,ripple_line,'Color','r')
if (export_as_pgf)
  plot_name = [devicename, '-phi_ripple_line'];
  PGFdata = [phi_ripple_line(:) - phi_line(1), ripple_line(:)];
  save([plot_name, plot_suf,'_pgf.dat'], 'PGFdata', '-ascii')
end
xlim([phi_start,phi_end])
xlabel('$\varphi_s$', 'Interpreter', 'latex')
ylabel('$1 / \hat B$', 'Interpreter', 'latex')

%% Observation points
% point = 1100;
%
% phi_l_point  = phi_bhat_line(point);
% bhat_point = bhat_line(point);
%
% hold on
% plot(phi_l_point, bhat_point, '.k', 'MarkerSize', 30)
% hold off
%
% iota    = mag_s.surface.tag_1.aiota;
% theta_0 = mag_s.fieldline.tag_1.xstart(3);
% phi_0   = mag_s.fieldline.tag_1.xstart(2);
%
% theta_point = theta_0 + iota * (phi_l_point - phi_0);
% phi_point   = phi_l_point;
%
% theta_point = mod(theta_point, 2*pi)
% phi_point   = phi_point%mod(phi_point, 2*pi)
