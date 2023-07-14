function [ds, sigma] = adaptSolveDiscs( ds, kernel, RHSfunc, opts )
% Adaptively solves Ksigma = RHS at discretization points on ds.
%
% IN: ds : a discs object
%       kernel : function handle, kernel being used in integration
%       RHSfunc : function handle used to compute RHS of BIE
% OUT : ds : discs object updated with new discretization points (if used)
%            sigma : array with solution density

if(nargin < 4)
    opts = [];
end

% Get the correct order if close region information provided in ds
% if not consider all panels as close


% zeroth refinement, initialize everything

Ninit = ds.chnkrs.npt; % initial number of discretization points
K = zeros(Ninit);





end