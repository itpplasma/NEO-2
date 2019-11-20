function b = b_fourier_cos(m, bm, theta)

  oldsize = size(theta);
  theta = theta(:).';

  b = sum(bm .* cos(m .* theta), 1);
  b = reshape(b, oldsize);

end
