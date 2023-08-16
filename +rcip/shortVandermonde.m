function V = shortVandermonde( x, m )
% *shortVandermonde* creates m columns of the Vandermonde matrix such that
% its columns are powers of the vector x (in increasing order)
%
% Syntax: V = shortVandermonde( x, m )
%
% Input:
%   x - initial vector
%   m - number of columns in the Vandermonde matrix to compute
%
% Output:
%   V - short Vandermonde matrix
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)
n = length(x);
x = reshape(x, n, 1); % x has to be a column vector
V = ones(n, m);
for i=2:m
    V(:, i) = V(:, i-1).*x;
end
end