function [r, d, d2] = pointsDisc(t, opts)

% Get points on disc
% return position, first, second derivatives, outward normals of 
% parameterized disc
%
% Inputs:
% t - points (in [0,2pi] to evaluate these quantities
%     (note that since we are dealing with identical discs, if we have 
%      nk discs in total then (j-1)*2pi/nk <= t <= j*2pi/nk for the
%      j-th disc)
%
% Optional inputs:
% R - float, radius of the disc (1)
% ctrs - float(2, ndiscs), x,y coordinates of center of disc ( [0,0] )

arguments
    t
    opts.R = 1.0
    opts.ctr = [0;0]
end

R = opts.R;
ctr = opts.ctr;

% Compute r, d, d2

ct = cos(t);
st = sin(t);

% fill in by hand the parametrization for the first circle

r = [R*ct + ctr(1); R*st + ctr(2)];
d = [-R*st; R*ct];
d2 = [-R*ct; -R*st];

end