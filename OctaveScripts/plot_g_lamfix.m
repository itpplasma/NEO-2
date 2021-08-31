function [ax, data, s] = plot_g_lamfix(h5file, quant, point, lam, x_cutoff, ax, style, color)

  if nargin < 8 || isempty(color) color='b';  end
  if nargin < 7 || isempty(style), style='-';  end
  if nargin < 6 || isempty(ax), ax = axes; end
  if nargin < 5 || isempty(x_cutoff), x_cutoff = 5; end

  epslam   = 1e-4;
  h5  = h52struct(h5file);

  if isa(point, 'char')
    points = cell(1);
    points{1} = point;
    colors = color;
  else
    points = point;
    colors = lines(numel(points));
  end
  num_p = numel(points);
  hold(ax, 'on')

  for k = 1:num_p
    h5p = h5.(points{k});
    g = interp1(h5p.lambda(:), h5p.(quant)(:, :).', lam, 'spline');

    % if any(abs(h5p.lambda - lam) < epslam)
    %   lami = find(abs(h5p.lambda(:) - lam) < epslam);
    %   lami = lami(1);
    %   h5p.lambda(lami)
    % else
    %     error('No lambda found!')
    % end

    Lx = h5p.x < x_cutoff;
    Lx(1) = false;

    %set(ax, 'YScale', 'log')
    plot(ax, h5p.x(Lx), g(Lx), 'LineWidth', 1.0, 'LineStyle', style, 'Color', colors(k,:));

    xlabel('x')
    ylabel('g')
    %xlabel('$\lambda$', 'Interpreter', 'Latex')
    %ylabel('$g_\mathrm{sp}$', 'Interpreter', 'Latex')
    %set(gca, 'XTick', [-1, -0.5, 0, 0.5, 1])
    %set(gca,'XMinorTick','on','YMinorTick','on')
    data.x = h5p.x(Lx);
    data.y = g(Lx);
  end
  hold(ax, 'off')
end
