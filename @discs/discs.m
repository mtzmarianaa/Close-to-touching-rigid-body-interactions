classdef discs
    % DISCS class which describes ndiscs non overlapping with same radius, 
    % different centers. They are represented by their values of its
    % position, first and second derivatives in parameters space
    % Since we also need to compute the arclength of those chunks we also
    % save the parametrizations
    %
    % Author: Mariana Martinez

    properties(Access = public)
        ctrs
        Rs
        listChnkrs
        chnkrs
        nDiscs
        nBreakPoints
    end

    methods
        function obj = discs(geom, p)
            % constructor for the discs class geom is the presets for the
            % geometry of the discs and p are the preferences for the
            % chunker objects
            if ( nargin < 1 )
                geom = [];
                geom.ctrs = [0;0];
                geom.Rs = 1;
                geom.nBreakPoints = 10;
            end
            if ( nargin < 2 )
                p = chunkerpref();
            else
                p = chunkerpref(p);
            end

            ctrs = geom.ctrs;
            nDiscs = size(ctrs, 2);
            obj.nDiscs = nDiscs;
            Rs = ones(1, nDiscs);
            nBreakPoints = 10*ones(1, nDiscs);

            if isfield(geom, 'Rs')
                Rs = geom.Rs;
            end
            if isfield(geom, 'nBreakPoints')
                if( size(geom.nBreakPoints, 2) < nDiscs)
                    nBreakPoints = geom.nBreakPoints(1)*ones(1, nDiscs);
                else
                    nBreakPoints = geom.nBreakPoints;
                end
            end

            obj.ctrs = ctrs;
            obj.Rs = Rs;
            obj.nBreakPoints = nBreakPoints;
            obj.nDiscs = nDiscs;
            listChnkrs(1, nDiscs) = chunker(); % list of nDiscs of chunkers
            for i=1:nDiscs
                breakpoints = linspace(0, 2*pi, nBreakPoints(i));
                center = ctrs(:, i);
                R = Rs(i);
                listChnkrs(i) = chunkerfuncbreakpoints( @(t) disc(t, center = center, radius = R), breakpoints, [], p );
            end
            chnkrs = merge(listChnkrs);
            obj.listChnkrs = listChnkrs;
            obj.chnkrs = chnkrs;
        end

    end

end