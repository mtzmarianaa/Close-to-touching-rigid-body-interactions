function [fvals] = f2_smooth(x)
% Test 2. Just returns sums of sines

fvals = 0.3*sin(x(1, :)) + 0.5*sin(x(2, :));

end