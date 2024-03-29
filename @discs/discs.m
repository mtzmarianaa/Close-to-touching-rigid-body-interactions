classdef discs
    % *discs* class which describes $N_{\Omega}=$nDiscs non overlaping
    % collection of discs with same radius but different centers. They are
    % represented by their values of their position, first and second
    % derivatives in parameter space. $\Gamma_k$ represent the
    % close-to-touching parts of the discs (if any) and $\Gamma_0$ is the
    % far region of the discs. Each close-to-toching region is represented
    % in both the coarse and fine mesh (necessary for the rcip method). 
    %
    % To do: automatization of process of finding the close-to-toching
    % region, make this more user friendly
    %
    % author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

    properties(Access = public)
        ctrs % center of the discs
        Rs % radii of discs
        listGammas % list chunker objects of gammas, close to touching parts
        listCoarseGammas % list chunker objects of gammas in the coarse discretization
        chnkrsGammas % chunker object with the information of all gammas
        nGammas % number of gammas, close to touching parts
        listFarChunks % list of chunkers (length of nDiscs) with information of far curve pieces
        gamma0 % chunkr object with information about the far part
        listChnkrs % list of chunkers with the information about the far and near curve pieces
        chnkrs % the chunker object with information from all discs
        nDiscs % number of discs
        nBreakPoints % number of breakpoints per disc
        infoClose % boolean, if information about close chunks is given or not
        I % if given, information about the chunks on each disc close to other discs
        I_closeReg % theta such that t +- theta defines the close region for each gamma in I
        indGammas % if close to touching information given, this is the index where gamma_i starts
        nB % limits for the K matrix
        nBCoarse % limits for the Kc matrix
    end

    properties(SetAccess = private)
        geom
        pClose
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
            ctrs = geom.ctrs;
            nDiscs = size(ctrs, 2);

            Rs = 0.75*ones(1, nDiscs); % all discs with same radii = 0.75 (NEVER USE 1)
            if isfield(geom, 'Rs') 
                Rs = geom.Rs;
                if( length(Rs) < nDiscs )
                    Rs = Rs*ones(1, nDiscs);
                end
            end
            geom.Rs = Rs;

            nBreakPoints = 10*ones(1, nDiscs);
            if isfield(geom, 'nBreakPoints')
                nBreakPoints = geom.nBreakPoints;
                if( size(geom.nBreakPoints, 2) < nDiscs)
                    nBreakPoints = geom.nBreakPoints(1)*ones(1, nDiscs);
                end
            end

            % Settings for the geometry of the discs
            obj.nDiscs = nDiscs;
            obj.geom = geom;
            obj.nDiscs = nDiscs;
            obj.ctrs = ctrs;
            obj.Rs = Rs;
            obj.nBreakPoints = nBreakPoints;

            if( nargin < 2)
                % Meaning that we dont have information about points close
                % to other discs, pClose is a structured array
                pClose = findpClose(obj);
                infoClose = true;
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

            obj.pClose = pClose;

            listGammas = [];
            listCoarseGammas = [];
            nGammas = 0;
            I = [];
            indGammas = 0;
            chnkrsGammas = [];
            listFarChunks = [];
            gamma0 = [];
            I_closeReg= [];
            nB = [];
            nBCoarse = [];

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
                        pClose(i).thetasReg = pi/6*ones(pClose(i).nClose, 1);
                    end

                    % Check that we have the minimum number of breakpoints
                    % on each disc
                    if( nBreakPoints(i) < 3*pClose(i).nClose )
                        nBreakPoints(i) = 3*pClose(i).nClose;
                    end

                end

                % Sort pClose, calculate amount of gammas and build I
                [I, nGammas, pClose] = buildI(obj);
                obj.I = I;
                obj.pClose = pClose; % Update

                if(isempty(I))
                    infoClose = false;
                end

            end

            if(infoClose)

                % Build list for gammas, close to touching interactiong and
                % the map for the neighbors
                [listGammas, neisMapClose, I_closeReg, listCoarseGammas] = buildGammas(obj);
                chnkrsGammas = merge(listGammas);

                % Build gamma0, neisMapFar, listFarChunks
                [gamma0, neisMapFar, listFarChunks] = buildGamma0(obj);

                % Add the missing neis
                indMissingClose = find(~chnkrsGammas.adj);
                indMissingFar = find(~gamma0.adj);
                % Merge
                chnkrsGammas = merge(listGammas);
                chnkrs = merge([gamma0, listGammas]);
                chnkrs .adj(2*gamma0.nch + indMissingClose) = neisMapFar;
                chnkrs.adj(indMissingFar) = gamma0.nch + neisMapClose;

                % Add the indices where the gammas start and nB
                indGammas = zeros(nGammas+2, 1);
                indGammas(1) = 0;
                indGammas(2) = gamma0.nch;
                nB = zeros(nGammas + 2, 1);
                nBCoarse = zeros(nGammas + 2, 1);
                nB(2) = gamma0.npt;
                nBCoarse(2) = gamma0.npt;
                for i=3:(nGammas+2)
                    indGammas(i) = indGammas(i-1) + listGammas(i-2).nch;
                    nB(i) = nB(i-1) + listGammas(i-2).npt;
                    nBCoarse(i) = nBCoarse(i-1) + listCoarseGammas(i-2).npt;
                end

            end

            
            if( ~ infoClose)
                % If we dont have close to touching parts
                nB = zeros( nDiscs + 1, 1);
                nB(1) = 0;
                for i=1:nDiscs
                    breakpoints = linspace(0, 2*pi, nBreakPoints(i));
                    center = ctrs(:, i);
                    R = Rs(i);
                    listChnkrs(i) = chunkerfuncbreakpoints( @(t) disc(t, center = center, radius = R), ...
                        breakpoints, [], p );
                    nB(i+1) = nB(i) + listChnkrs(i).npt;
                end
                chnkrs = merge(listChnkrs);
            end

            % Add these properties to the object
            obj.infoClose = infoClose;
            obj.chnkrsGammas = chnkrsGammas;
            obj.chnkrs = chnkrs;
            obj.listGammas = listGammas;
            obj.listCoarseGammas = listCoarseGammas;
            obj.nGammas = nGammas;
            obj.listFarChunks = listFarChunks;
            obj.gamma0 = gamma0;
            obj.I = I;
            obj.indGammas = indGammas;
            obj.nB = nB;
            obj.nBCoarse = nBCoarse;
            obj.I_closeReg = I_closeReg;

            if( infoClose )
                listChnkrs = [gamma0, listGammas];
            end

            obj.listChnkrs = listChnkrs;
                         

        end

    end

end


