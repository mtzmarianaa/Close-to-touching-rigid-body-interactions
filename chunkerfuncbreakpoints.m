function chnkr = chunkerfuncbreakpoints(fcurve, ts, cparams, pref)

if ( nargin < 3 )
    cparams = [];
end

if ( nargin < 4 )
    pref = chunkerpref();
else
    pref = chunkerpref(pref);
end

ta = 0; tb = 2*pi;
if ( isfield(cparams,'ta') ), ta = cparams.ta; end
if ( isfield(cparams,'tb') ), tb = cparams.tb; end

if ( any(ts < ta | ts > tb) )
    error('Panel breakpoints must lie within parameter domain.');
end

nch = length(ts)-1;
k = pref.k;
x = lege.exps(k);
dim = checkcurveparam(fcurve, ta);
pref.dim = dim;

chnkr = chunker(pref);
chnkr = chnkr.addchunk(nch);
for i = 1:nch
    a = ts(i);
    b = ts(i+1);
    tt = a + (b-a)*(x+1)/2;
    [z, dz, dzz] = fcurve(tt);
    chnkr.r(:,:,i)  = reshape(z,   dim, k);
    chnkr.d(:,:,i)  = reshape(dz,  dim, k);
    chnkr.d2(:,:,i) = reshape(dzz, dim, k);
    chnkr.h(i) = (b-a)/2;
    chnkr.adj(:,i) = [i-1 i+1];
end
chnkr.adj(1,1)   = nch;
chnkr.adj(2,nch) = 1;
chnkr.n = normals(chnkr);

end
