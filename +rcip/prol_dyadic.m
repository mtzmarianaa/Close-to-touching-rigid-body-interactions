function P = prol_dyadic( x, nRef )
% Given x ( n Gauss-Legendre nodes on [-1, 1]) we build the prolongation
% matrix that interpolates data from this discretization to a finer one. 
% This finer discretization was done using nRef dyadic refinements of the
% original coarse discretization. Basically solving for P with two
% Vandermonde matrices.

n = length(x); % Get the number of Gauss-Legendre nodes on the coarse discretization

% Get the Legendre nodes in the finer discretization
n_fine = 2*(nRef)*n;
x_fine = zeros( n_fine, 1 );
x_fine( (nRef + 1)*n + 1, 1) = (x(1)+x(end))/2;
% Fill in x_fine
for i=1:(nRef-1)
    a = 2^(-i);
    b = 2^(-i+1);
    shifted = a + (b-a)*(x+1)/2;
    x_fine( ((i-1)*n + 1 ):i*n ) = - flip(shifted);
    x_fine( (n_fine-(i)*n + 1):( n_fine-(i-1)*n ) ) = shifted;
end

% Fill the middle part
a = 0;
b = 2^(-nRef + 1);
shifted = a + (b-a)*(x+1)/2;
x_fine( ((nRef-1 )*n + 1):((nRef +1 )*n) ) = [ -flip(shifted); shifted ];


% Build the Vandermonde matrices
C = rcip.shortVandermonde(x, n);
F = rcip.shortVandermonde(x_fine, n);

P = F/C;



end