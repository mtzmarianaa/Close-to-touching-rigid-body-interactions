function submat = prol2( x )
% Given x ( n Gauss-Legendre nodes on [-1, 1]) we build the prolongation
% matrix that interpolates data from this discretization to a finer one
% with 2*n Gauss-Legendre nodes. Basically solving for P with two
% Vandermonde matrices.

n = length(x); % Get number of Gauss-Legendre nodes on the coarse discretization


% Build the Vandermonde matrices
x_fine = [x-1; x+1]; % Legendre nodes in the finer discretization
x_fine = x_fine./2;
C = shortVandermonde(x, n);
F = shortVandermonde(x_fine, n);

submat = F/C;

end

