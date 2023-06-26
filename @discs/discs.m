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
        pClose
        infoClose
    end

    methods
        function obj = discs(geom, pClose, p)
            % constructor for the discs class geom is the presets for the
            % geometry of the discs and p are the preferences for the
            % chunker objects
            if ( nargin < 1 )
                geom = [];
                geom.ctrs = [0;0];
                geom.Rs = 1;
                geom.nBreakPoints = 10;
            end
            if( nargin < 2)
                % Meaning that we dont have information about points close
                % to other discs, pClose is a structured array
                pClose = [];
                infoClose = false;
            else
                % See if we have the correct amount of discs 
                if( length(pClose) ~= size(geom.ctrs, 2) )
                    pClose = [];
                    infoClose = false;
                else
                    infoClose = true;
                end
            end
            if ( nargin < 3 )
                p = chunkerpref();
            else
                p = chunkerpref(p);
            end

            % Settings for the geometry of the discs
            ctrs = geom.ctrs;
            nDiscs = size(ctrs, 2);
            obj.nDiscs = nDiscs;
            Rs = 0.75*ones(1, nDiscs); % all discs with same radii = 0.75 (NEVER USE 1)
            nBreakPoints = 10*ones(1, nDiscs);

            if isfield(geom, 'Rs')
                % No field for radii, standard option
                Rs = geom.Rs;
            end
            if isfield(geom, 'nBreakPoints')
                % No field for number of breakpoints, standard option
                nBreakPoints = geom.nBreakPoints;
                if( size(geom.nBreakPoints, 2) < nDiscs)
                    % If number of breakpoints is not correct, set all discs
                    % with the same radius
                    nBreakPoints = geom.nBreakPoints(1)*ones(1, nDiscs);
                end
            end

            % Settings for the geometry of the close to touching region of
            % the discs
            if( infoClose  )
                % Meaning that we DO have close points

                for i=1:nDiscs
                    % Check that we have the minimum number of breakpoints
                    % on each disc
                    if( nBreakPoints(i) < 3*pClose(i).nClose )
                        nBreakPoints(i) = 3*pClose(i).nClose;
                    end

                    % If we aren't given the number of close points,
                    % calculate them
                    if ~isfield(pClose(i), 'nClose')
                        pClose.nClose = length(pClose(i).thetas);
                    end

                    % If we aren't given the angle for the region to be
                    % considered as close, set it to pi/3
                    if ( ~isfield(pClose(i), 'thetasReg') || length(pClose(i).thetasReg) ~=  pClose(i).nClose)
                        pClose(i).thetasReg = pi/6*ones(pClose(i).nClose);
                    end

                    % Sort them
                    [pClose(i).thetas, I] = sort(pClose(i).thetas);
                    pClose(i).thetas = mod(pClose(i).thetas, 2*pi); % Just that everything agrees
                    pClose(i).thetasReg = pClose(i).thetasReg(I);

                end
            end

            % Set the parameters we already have
            obj.ctrs = ctrs;
            obj.Rs = Rs;
            obj.nBreakPoints = nBreakPoints;
            obj.nDiscs = nDiscs;
            obj.pClose = pClose;
            obj.infoClose = infoClose;
            listChnkrs(1, nDiscs) = chunker(); % list of nDiscs of chunkers

            % For each disc we need to add the points that are close AND
            % pointClose +- pi/6. Then add the other breakpoints
            % accordingly. Remember that pClose.thetasClose is in the angle
            % space, not in the coordinate space (so this is easier)

            if( infoClose )
                % If we have close to touching points
                for i=1:nDiscs
                    center = ctrs(:, i);
                    R = Rs(i);
                    thisnClose = pClose(i).nClose;
                    breakpoints = zeros(3*thisnClose, 1); % breakpoints to initite chunkie
                    farIntervals = zeros(thisnClose, 2); % intervals not in the close region
                    lenFarIntervals = zeros(thisnClose, 1); 
                    % Add the close to touching points and their regions
                    for k=1:thisnClose
                        j = 3*k;
                        thetaR = pClose(i).thetasReg(k);
                        thetaC = pClose(i).thetas(k);
                        % add the close region
                        breakpoints(j-2) = thetaC - thetaR/2;
                        breakpoints(j-1) = thetaC;
                        breakpoints(j) = thetaC + thetaR/2;
                    end
                    breakpoints = mod(breakpoints, 2*pi);
                    % Compute the far intervals
                    for k=1:(thisnClose - 1)
                        j = 3*k;
                        farIntervals(k, :) = [breakpoints(j-2), breakpoints(j)]; 
                        lenFarIntervals(k) = breakpoints(j) - breakpoints(j-2);
                    end
                    % For the last one
                    farIntervals(thisnClose, :) = [breakpoints(3*thisnClose),  breakpoints(1) ]; 
                    lenFarIntervals(thisnClose) = breakpoints(1) - breakpoints(3*thisnClose);
                    farIntervals = mod(farIntervals, 2*pi);
                    lenFarIntervals = mod(lenFarIntervals, 2*pi);
                    
                    
                    % Then determine the length of all the far section and
                    % amount of breakpoints to put here
                    lenFarSection = sum(lenFarIntervals);
                    nBreakPointsFar = round( (nBreakPoints(i) - 3*thisnClose )*lenFarIntervals./lenFarSection );
                    breakpointsFar = zeros( sum(nBreakPointsFar), 1 );
                    nBk = [0; nBreakPointsFar];
                    for k=2:(thisnClose + 1)
                        nBk(k) = nBk(k - 1) + nBk(k);
                    end
                    
                    % Compute the actual breakpoints
                    for k = 1:thisnClose
                        if( farIntervals(k, 1) < farIntervals(k, 2) )
                            grid = linspace( farIntervals(k, 1), farIntervals(k, 2), nBreakPointsFar(k) + 2  )';
                        else
                            grid = linspace( farIntervals(k, 1), 2*pi+farIntervals(k, 2), nBreakPointsFar(k) + 2  )';
                            grid = mod(grid, 2*pi);
                        end
                        breakpointsFar( (nBk(k)+1):nBk(k+1) ) = grid(2:(length(grid)-1));
                    end

                    % Add them, sort them, put them in the interval 0 2pi
                    breakpoints = [breakpoints; breakpointsFar];
                    breakpoints = sort( mod(breakpoints, 2*pi) );
                    % Make sure everything is connected
                    breakpoints = [0; 2*pi; breakpoints];
                    breakpoints = unique(breakpoints);

                    % Add this info to the list of chunkers
                    listChnkrs(i) = chunkerfuncbreakpoints( @(t) disc(t, center = center, radius = R), breakpoints, [], p );

                end
            else
                % If we dont have close to touching parts
                for i=1:nDiscs
                    breakpoints = linspace(0, 2*pi, nBreakPoints(i));
                    center = ctrs(:, i);
                    R = Rs(i);
                    listChnkrs(i) = chunkerfuncbreakpoints( @(t) disc(t, center = center, radius = R), breakpoints, [], p );
                end
            end
            % Merge into a single chunkr object
            chnkrs = merge(listChnkrs);
            obj.listChnkrs = listChnkrs;
            obj.chnkrs = chnkrs;

        end

    end

end

