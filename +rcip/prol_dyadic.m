function P = prol_dyadic( k, nRef )
% Given 4 panels with k Legendre nodes each we build the prolongation
% matrix that maps those nodes to nodes on the finer discretization
% (dyadically refined towards the middle)
% Matrix is of the form size 4k x ( 2k + 2k(nRef + 1) ):
%    I   0   0   0
%    0   P1  0   0
%    0   0   P2  0
%    0   0    0   I


% Use the functions defined below

P1 = pol_dyadicPart1(k, nRef);
P2 = pol_dyadicPart2(k, nRef);
Ik = speye(k);

P = blkdiag( Ik, P1, P2, Ik );

end


function P1 = pol_dyadicPart1( k, nRef )
% Part of the prolongation matrix where the dyadic refinement is towards 1

x_lege = lege.exps(k);

n_fine = k*(nRef + 1);
x_fine = zeros(n_fine, 1);

% Fill the first part
a = -1;
if( nRef == 0 )
    b = 1;
else
    b = 0;
end
shifted = a + (b-a)*(x_lege+1)/2;
x_fine( 1:k  ) = shifted;

breakPoints = ones(nRef, 1);
for i=1:(nRef-1)
    breakPoints(i) = 1 - 2^(-i);
end

% Fill x_fine
for i=1:nRef
    a = b; % Start where we ended up
    b = breakPoints(i);
    shifted = a + (b-a)*(x_lege+1)/2;
    x_fine( (i*k + 1):(i+1)*k  ) = shifted;
end

% Build the Vandermonde matrices
C1 = rcip.shortVandermonde(x_lege, k);
F1 = rcip.shortVandermonde(x_fine, k);

P1 = F1/C1;

end


function P2 = pol_dyadicPart2( k, nRef )
% Part of the prolongation matrix where the dyadic refinement is towards -1

x_lege = lege.exps(k);

n_fine = k*(nRef + 1);
x_fine = zeros(n_fine, 1);

% Fill the first part
b = 1;
if( nRef == 0 )
    a = -1;
else
    a = 0;
end
shifted = a + (b-a)*(x_lege+1)/2;
x_fine( (k*nRef + 1):( k*(nRef + 1) )  ) = shifted;

breakPoints = -ones(nRef, 1);
for i=1:(nRef-1)
    breakPoints(i) = -1 + 2^(-i);
end

% Fill x_fine
for i=1:nRef
    b = a; % End where we started last time
    a = breakPoints(i); % Start before
    shifted = a + (b-a)*(x_lege+1)/2;
    x_fine( (k*(nRef - i) + 1):(k*(nRef - i + 1))  ) = shifted;
end

% Build the Vandermonde matrices
C2 = rcip.shortVandermonde(x_lege, k);
F2 = rcip.shortVandermonde(x_fine, k);

P2 = F2/C2;

end


