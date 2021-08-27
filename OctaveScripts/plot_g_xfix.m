function [ax, data, s] = plot_g_xfix(h5file, quant, point, x, ax, style, color)

  if nargin < 7 || isempty(color), color='b';  end
  if nargin < 6 || isempty(style), style='-';  end
  if nargin < 5 || isempty(ax), ax = axes; end

  epsx   = 1e-6;

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
    h5p.lambda(:);

    if any(abs(h5p.x - x) < epsx)
      kx = find(abs(h5p.x(:) - x) < epsx);
      kx = kx(1);
    else
      error('No velocity found!')
    end

    plot(ax, h5p.lambda, h5p.(quant)(kx,:), 'LineWidth', 1.0, 'LineStyle', style, 'Color', colors(k,:));
    xlabel('\lambda')
    ylabel('g')
    %xlabel('$\lambda$', 'Interpreter', 'Latex')
    %ylabel('$g_\mathrm{sp}$', 'Interpreter', 'Latex')
    %set(gca, 'XTick', [-1, -0.5, 0, 0.5, 1])
    %set(gca,'XMinorTick','on','YMinorTick','on')
    data.x =  h5p.lambda;
    data.y = h5p.(quant)(kx,:);
    hold(ax, 'off')
  end
end
