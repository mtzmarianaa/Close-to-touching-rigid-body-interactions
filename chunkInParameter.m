function ch = chunkInParameter(chnkr, tha, thb)
% Given a chunker object and the ends points in parameter space
% get that chunk in the parameter space
%
% Input:
%        chnkr - the chunker object
%        tha  - point at the beggining of the chunk in parameter space
%        thb  - point at the end of the chunk in parameter space

x = lege.exps(chnkr.k);
ch = tha + (thb - tha)*(x + 1)/2;

end






