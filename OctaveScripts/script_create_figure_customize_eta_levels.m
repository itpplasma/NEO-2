N = 100;
x = 0:N-1;
n = x.^2;
alpha = 0.5;
t = tanh(2*alpha*n./(N-1)^2)*(N-1)^2./tanh(2*alpha);

c_temp = cos(pi*n./(N-1)^2 - pi);
b = sign(c_temp) .* abs(c_temp).^alpha;
c = (b - min(b)) *(N-1)^2 ./(max(b) - min(b));

plot(x, n, 'kx-', x, t, 'b+-', x, c, 'r*-');
xlabel('index');
ylabel('eta');
legend('1 identity', '2 tanh', '3 cos');
title(['\alpha = ', num2str(alpha)]);

print(['documentation_customize_eta_levels.eps'], '-depsc');
print(['documentation_customize_eta_levels.png'], '-dpng');
