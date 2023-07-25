function [R, G] = buildR(chunkerCoarse, chunkerFine, K22_inv, P,  kernel)
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

nRef = floor( chunkerFine.nch/4 - 2  );

if(nargin < 4)
    P = rcip.prol_dyadic(chunkerCoarse.k, nRef);
    P = blkdiag(P, P);
end

if(nargin < 3)
    K22 = chunkermat(chunkerFine, kernel);
    K22_inv = inv(K22);
end

% Build R and G if needed
 if( nargout < 2) 
     R = Wc*P'*Wf*K22_inv*P;
 else
     G = Wc*P'*Wf;
     R = G*K22_inv*P;
 end


end