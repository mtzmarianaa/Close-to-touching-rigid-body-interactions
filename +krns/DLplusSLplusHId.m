function K = DLplusSLplusHId(s,t)

K = krns.DL_kern(s,t) + krns.SL_kern(s,t);

if( isequal( s.r, t.r) )
    K = K + 0.5*eye(size(K));
end

end