%helper function that performs decoy state analysis
function [Y1L,Y1U] = decoyIndependent(expecJointTest,decoys,probListAlice,probsDecoy,ptest,t,mu,options)
    arguments
        expecJointTest (:,:,:) double {mustBeGreaterThanOrEqual(expecJointTest,0)};
        decoys (:,:) cell;        
        probListAlice (:,1) double {mustBeProbDist}
        probsDecoy (:,1) double {mustBeProbDist}
        ptest (1,1) double {mustBeInRange(ptest,0,1),mustBeSameTotalProb(ptest,expecJointTest)}
        t (:,:,:) double {mustBeGreaterThanOrEqual(t,0),SameNumberOfDecoysAndPages(decoys,t)}
        mu (:,:,:) double {mustBeGreaterThanOrEqual(mu,0),SameNumberOfDecoysAndPages(decoys,mu)}
        options.decoyTolerance (1,1) double {mustBeGreaterThanOrEqual(options.decoyTolerance,0)}= 1e-12;
        options.ChoiTolerance (1,1) double {mustBeGreaterThanOrEqual(options.ChoiTolerance,0)}= 1e-10;
        options.decoySolver (1,1) string = "sdpt3";
        options.decoyPrecision (1,1) string = "high";
        options.photonCutOff (1,1) double {mustBeInteger,mustBePositive} = 10;
        options.verboseLevel (1,1) double {mustBeInteger,mustBeNonnegative} = 1;
    end
   
    %number of intensities used
    n_decoy=numel(decoys);

    %m = #number of rows of observations
    %n = #number of colums of observations
    m = size(expecJointTest,1);
    n = size(expecJointTest,2);
    
    %Empty vectors for upper and lower bounds
    Y1L = zeros(m,n);
    Y1U = ones(m,n);
    
    %cut-off for decoy
    n_photon= options.photonCutOff;
    
    %Number of intensities
    n_signal_int = size(decoys{1},2);

    %generate matrix of decoy intensities
    decoyMat = cell2mat(decoys);

    % %Possonian distribution
    probDist = @(intensity,n) exp(-intensity).*intensity.^n./factorial(n);

    %Thermal distribution
    % probDist=@(meanPh,n) meanPh.^n/(1+meanPh).^(1+n);
    
    %Tolerance in Decoy SDP
    decoy_tolerance = options.decoyTolerance;
    choi_tolerance = options.ChoiTolerance;

    %Pick n=0,1,...-photon components
    M = {1,n_photon+1};
    for i = 1:n_photon+1
        Mi = zeros(m, m*(n_photon + 1));
        indr = [1:m]';
        indc = [i:n_photon+1:m*(n_photon + 1)]';
        indx = sub2ind(size(Mi),indr,indc);
        Mi(indx) = 1;
        M{i}=Mi;
    end

    %Calculate upper and lower bounds on the expectations
    upperbndsExpec = ones(m,n,numel(decoys));
    lowerbndsExpec = zeros(m,n,numel(decoys));

    for kdecoy = 1:numel(decoys)
        for row = 1:numel(probListAlice)
            upperbndsExpec(row,:,kdecoy) = min(max((expecJointTest(row,:,kdecoy) + t(row,:,kdecoy) + mu(row,:,kdecoy))/(probListAlice(row)*probsDecoy(kdecoy)*ptest),0),1);
            lowerbndsExpec(row,:,kdecoy) = min(max((expecJointTest(row,:,kdecoy) - t(row,:,kdecoy) - mu(row,:,kdecoy))/(probListAlice(row)*probsDecoy(kdecoy)*ptest),0),1);
        end
    end
    
    
    %solve for upper bound 1-photon
    for i=1:m
        for j=1:n
            try
                %select correct intensities
                decoyList = decoyMat(:,i);

                %Set up SDP for upper bounds
                % Y is the matrix representing all yields
                % J0 is the Choi matrix for the 0-photon component
                % J1 is the Choi matrix for the 1-photon component
                % In both Choi matrices Bob's system is first and Alice's
                % system is second to satisfy e.g. Y_0 = Tr[J_0 (F_y otimes (rho_x)^T)] 
                cvx_begin quiet
                    %CVX solver
                    cvx_solver(convertStringsToChars(options.decoySolver));
                    variable Y(n_photon+1,1)
                    maximize Y(2)
                    subject to
                        % 0 <= Y <=1
                        vec(Y) >= vec(zeros(n_photon+1,1));
                        vec(Y) <= vec(ones(n_photon+1,1));             
                        
                        %Usual decoy bounds rewritten as matrix times vector
                        for kdecoy = 1:n_decoy
                            %total probability for given intensity to have
                            %<= n_photon photons
                            %create empty vector which is filled for each
                            %intensity again

                            %select subintensity
                            intensity = decoyList(kdecoy);
                            %calculate prob. distribution
                            Pmu = probDist(intensity,0:1:n_photon);
                            %create list with all subintensities
                            pmu_tot = sum(Pmu);
                            
                            %Decoy constraints
                            Pmu*Y >= lowerbndsExpec(i,j,kdecoy) - (1-pmu_tot) - decoy_tolerance;
                            Pmu*Y <= upperbndsExpec(i,j,kdecoy) + decoy_tolerance;                   
                        end                        
                cvx_end

            if options.verboseLevel>=2
                disp(cvx_status)
            end

            if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end
            catch error 
              rethrow(error);%  fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end

            %Store upper bound and make sure 0 <= Y1U <= 1 
            % (can be violated due to numerics) 
            Y1U(i,j) = min(max(Y(2),0),1);            
        end
    end
        
    %solve for lower bounds
    %solve for upper bound 1-photon
    for i=1:m
        for j=1:n
            try
                %select correct intensities
                decoyList = decoyMat(:,i);

                %Set up SDP for upper bounds
                % Y is the matrix representing all yields
                % J0 is the Choi matrix for the 0-photon component
                % J1 is the Choi matrix for the 1-photon component
                % In both Choi matrices Bob's system is first and Alice's
                % system is second to satisfy e.g. Y_0 = Tr[J_0 (F_y otimes (rho_x)^T)] 
                cvx_begin quiet
                    %CVX solver
                    cvx_solver(convertStringsToChars(options.decoySolver));
                    variable Y(n_photon+1,1)
                    minimize Y(2)
                    subject to
                        % 0 <= Y <=1
                        vec(Y) >= vec(zeros(n_photon+1,1));
                        vec(Y) <= vec(ones(n_photon+1,1));             
                        
                        %Usual decoy bounds rewritten as matrix times vector
                        for kdecoy = 1:n_decoy
                            %total probability for given intensity to have
                            %<= n_photon photons
                            %create empty vector which is filled for each
                            %intensity again

                            %select subintensity
                            intensity = decoyList(kdecoy);
                            %calculate prob. distribution
                            Pmu = probDist(intensity,0:1:n_photon);
                            %create list with all subintensities
                            pmu_tot = sum(Pmu);
                            
                            %Decoy constraints
                            Pmu*Y >= lowerbndsExpec(i,j,kdecoy) - (1-pmu_tot) - decoy_tolerance;
                            Pmu*Y <= upperbndsExpec(i,j,kdecoy) + decoy_tolerance;                   
                        end                        
                cvx_end

            if options.verboseLevel>=2
                disp(cvx_status)
            end

            if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end
            catch error 
              rethrow(error);%  fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end

            %Store upper bound and make sure 0 <= Y1U <= 1 
            % (can be violated due to numerics) 
            Y1L(i,j) = min(max(Y(2),0),1);          
        end
    end
end

%% validation functions
function SameNumberOfDecoysAndPages(decoys,conditionalExpectations)
if numel(decoys) ~= size(conditionalExpectations,3)
    throwAsCaller(MException("decoyAnalysisIndependentLP:DecoysAndPagesDontMatch",...
        "The number of decoy intensities does not match the number of pages for conditionalExpectations."))
end
end

function mustBeSameTotalProb(totProb,dist)
if ~ismembertol(sum(dist,"all"),totProb)
    throw(MException("decoyAnalysisConstrFinite:TotalTestProbViolated",...
        "The joint distribution sum to the probability of testing up to tolerance."))
end
end