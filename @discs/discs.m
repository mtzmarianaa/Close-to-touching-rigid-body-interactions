classdef discs
    % DISCS class which describes ndiscs non overlapping with same radius, 
    % different centers. They are represented by their values of its
    % position, first and second derivatives in parameters space
    % Since we also need to compute the arclength of those chunks we also
    % save the parametrizations
    %
    % Author: Mariana Martinez

    properties(Access = public)
        ctrs % center of the discs
        Rs % radii of discs
        listGammas % list of gammas, close to touching parts
        chnkrsGammas % chunker object with the information of all gammas
        nGammas % number of gammas, close to touching parts
        listFarChunks % list of chunkers (length of nDiscs) with information of far curve pieces
        gamma0 % chunkr object with information about the far part
        chnkrs % the chunker object with information from all discs
        nDiscs % number of discs
        nBreakPoints % number of breakpoints per disc
        infoClose % boolean, if information about close chunks is given or not
        I % if given, information about the chunks on each disc close to other discs
        indGamma % if close to touching information given, this is the index where the close to touching part
                         % starts in chnkrs (i.e. where K22 begins in K)
    end

    methods
        function obj = discs(geom, pClose, p)
            % constructor for the discs class geom is the presets for the
            % geometry of the discs and p are the preferences for the
            % chunker objects.
            %
            % Input:
            %        geom - geometry options for the discs (at least ctrs provided)
            %                   geom.ctrs = centers of the discs
            %                   geom.Rs = radii of the discs
            %                   geom.nBreakPoint = number of breakpoints to
            %                           use on each disc, if just one, all discs
            %                           with same number of breakpoints
            %        pClose - information about the close to touching
            %            regions (at least thetas, discClose, pRef given)
            %                     pClose(i).data = (nClose x 3 matrix) at disc i the point
            %                     closest to another disc in parameter
            %                     space. The first column is the
            %                     point in parameter (angle) space. The 
            %                     second column is the index of
            %                     the other disc which that point is
            %                     closest to. The third column is the index
            %                     of the point in the other disc that is
            %                     closest to disc i.
            %                     pClose(i).nClose = number of points at
            %                     disc i considered to be in the close
            %                     region
            %                     pClose(i).thetasReg = angle of region
            %                     considered to be close
            %        infoClose - boolean, indicates if we have information
            %                         regarging close to touching
            %                         interactions.
            %        p - chunker preferences
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
            thetas = []; % Here is where we are going to save (or not) the parametrization thetas

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

            listGammas = [];
            nGammas = 0;
            I = [];
            indGamma = 0;
            chnkrsGammas = [];

            % Settings for the geometry of the close to touching region of
            % the discs
            if( infoClose  )
                % Meaning that we DO have close points

                for i=1:nDiscs

                    % If we aren't given the number of close points,
                    % calculate them
                    if (~isfield(pClose(i), 'nClose') || isempty(pClose(i).nClose) )
                        pClose(i).nClose = length(pClose(i).thetas);
                    end

                    % If we aren't given the angle for the region to be
                    % considered as close, set it to pi/3
                    if ( ~isfield(pClose(i), 'thetasReg') || ...
                            length(pClose(i).thetasReg) ~=  pClose(i).nClose)
                        pClose(i).thetasReg = pi/3*ones(pClose(i).nClose, 1);
                    end

                    % Check that we have the minimum number of breakpoints
                    % on each disc
                    if( nBreakPoints(i) < 3*pClose(i).nClose )
                        nBreakPoints(i) = 3*pClose(i).nClose;
                    end

                end

                % Sort pClose, calculate amount of gammas and build I
                [I, nGammas, pClose] = dsc.buildI(pClose, nDiscs);

                % Build list for gammas, close to touching interactiong and
                % the map for the neighbors
                [listGammas, neisMapClose] = dsc.buildGammas(geom, pClose, I, nDiscs);
                chnkrsGammas = merge(listGammas);

                % Build gamma0, neisMapFar, listFarChunks
                [gamma0, neisMapFar, listFarChunks] = dsc.buildGamma0(geom, pClose, I, nDiscs);

                % Add the missing neis
                indMissingClose = find(~chnkrsGammas.adj);
                indMissingFar = find(~gamma0.adj);
                % Merge
                chnkrs = merge([gamma0, listGammas]);
                chnkrs .adj(2*gamma0.nch + indMissingClose) = neisMapFar;
                chnkrs.adj(indMissingFar) = gamma0.nch + neisMapClose;
                
                indGamma = gamma0.nch + 1;

            end

            
            if( ~ infoClose)
                % If we dont have close to touching parts
                for i=1:nDiscs
                    breakpoints = linspace(0, 2*pi, nBreakPoints(i));
                    center = ctrs(:, i);
                    R = Rs(i);
                    listChnkrs(i) = chunkerfuncbreakpoints( @(t) disc(t, center = center, radius = R), ...
                        breakpoints, [], p );
                end
                chnkrs = merge(listChnkrs);
            end

            % Add these properties to the object
            obj.chnkrs = chnkrs;
            obj.listGammas = listGammas;
            obj.nGammas = nGammas;
            obj.listFarChunks = listFarChunks;
            obj.gamma0 = gamma0;
            obj.I = I;
            obj.indGamma = indGamma;
                         

        end

    end

end


