function b = b_fourier_sin(m, bm, theta)

    oldsize = size(theta);
    theta = theta(:).';

    b = sum(bm .* sin(m .* theta), 1);
    b = reshape(b, oldsize);

end
