function [xknots] = compute_knots(nrknots, phi_x_max, bspline_dist)
  gam_all = 0.0;
  for k = 1:nrknots-1
    gam_all = gam_all + bspline_dist.^k;
  end

  xknots(1) = 0.0;
  x_del = phi_x_max / gam_all;
  for k = 1:nrknots-1
    xknots(k+1) = xknots(k) + x_del * bspline_dist**k;
  end
end
