nrknots = 5;
phi_x_max = 4;

xknots05 = compute_knots(nrknots, phi_x_max, 0.5);
xknots10 = compute_knots(nrknots, phi_x_max, 1.0);
xknots15 = compute_knots(nrknots, phi_x_max, 1.5);
xknots15b = compute_knots(2*nrknots-1, phi_x_max, 1.5);
xknots20 = compute_knots(nrknots, phi_x_max, 2.0);
xknots20b = compute_knots(2*nrknots-1, phi_x_max, 2.0);

plot(xknots05, zeros(size(xknots05)), 'k*',...
     xknots10, 1*ones(size(xknots10)), 'k+',...
     xknots15, 2*ones(size(xknots15)), 'kx',...
     xknots15b, 2.5*ones(size(xknots15b)), 'kv',...
     xknots20, 3*ones(size(xknots20)), 'ko',...
     xknots20b, 3.5*ones(size(xknots20b)), 'k^');
xlabel('eta')
legend('bspline_dist = 0.5', ' = 1.0', ' = 1.5', ' = 1.5', ' = 2.0', ' = 2.0');

print(['documentation_bspline_distance.eps'], '-depsc');
print(['documentation_bspline_distance.png'], '-dpng');
