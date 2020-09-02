% \brief Make input(array) even.
%
% At the moment this is done by 'rounding' towards zero (fix).
function output = make_even(in)
  output = 2*fix(in/2);
end
