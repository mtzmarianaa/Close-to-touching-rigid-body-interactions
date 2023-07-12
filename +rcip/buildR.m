function R = buildR(chunkerCoarse, chunkerFine, kernel, K22_inv)
% Build the matrix    R = (W−1c )*( P' )*( Wf )*(  K−122 )*( P )

% We need Wc, Wf, P, K22
Wc = weights(chunkerCoarse);
Wc = 1./(Wc(:));
Wc = spdiags(Wc);
nc = chunkerCoarse.npt;
Wc = speye(nc).*Wc;

Wf = weights(chunkerFine);
Wf = Wf(:);
Wf = spdiags(Wf);
nf = chunkerFine.npt;
Wf = speye(nf).*Wf;

nRef = floor( chunkerFine.nch/2 );
x = lege.exps( chunkerCoarse.k );
P = rcip.prol_dyadic(x, nRef);
P = [P P];

if(nargin < 4)
    K22 = chunkermat(chunkerFine, kernel);
    K22_inv = inv(K22);
end

% Build R
R = Wc*P'*Wf*K22_inv*P;


end