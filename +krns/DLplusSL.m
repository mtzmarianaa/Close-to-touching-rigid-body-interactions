function K = DLplusSL(s,t)

K = krns.DL_kern(s,t) + krns.SL_kern(s,t);


end