phi_b = m0*T-n0*P;
phi_b_per = mod(phi_b, 2*pi);
s_surf = interp1(phi_vpt, s_resc(ir,:), phi_b_per, 'linear', 'extrap');
theta_surf = interp1(phi_vpt, phi_resc(ir,:)-phi_vpt, phi_b_per, 'linear', 'extrap')/m0+T;
