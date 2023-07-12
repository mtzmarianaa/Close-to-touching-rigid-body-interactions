function V = shortVandermonde( x, m )
% Creates n columns of the Vandermonde matrix such that its columns are 
% powers of the vector v (in increasing order)

n = length(x);
x = reshape(x, n, 1); % x has to be a column vector
V = ones(n, m);

for i=2:m
    V(:, i) = V(:, i-1).*x;
end

end