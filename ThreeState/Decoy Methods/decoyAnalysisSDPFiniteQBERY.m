%helper function that performs decoy state analysis
function [Y1L,Y1U] = decoyAnalysisSDPFiniteQBERY(expecJointTest,statesAlice,POVMsBob,decoys,probListAlice,probsDecoy,ptest,photonBound,t,mu,options)
    arguments
        expecJointTest (:,:,:) double {mustBeGreaterThanOrEqual(expecJointTest,0)}
        statesAlice (1,:) cell;
        POVMsBob (1,:) cell;
        decoys (:,1) double {mustBeGreaterThanOrEqual(decoys,0),SameNumberOfDecoysAndPages(decoys,expecJointTest)}        
        probListAlice (:,1) double {mustBeProbDist}
        probsDecoy (:,1) double {mustBeProbDist}
        ptest (1,1) double {mustBeInRange(ptest,0,1),mustBeSameTotalProb(ptest,expecJointTest)}
        photonBound (:,2) double {mustBeNonnegative,mustBeLessThanOrEqual(photonBound,1)}
        t (:,:,:) double {mustBeGreaterThanOrEqual(t,0),SameNumberOfDecoysAndPages(decoys,t)}
        mu (:,:,:) double {mustBeGreaterThanOrEqual(mu,0),SameNumberOfDecoysAndPages(decoys,mu)}
        options.decoyTolerance (1,1) double {mustBeGreaterThanOrEqual(options.decoyTolerance,0)}= 1e-14;
        options.ChoiTolerance (1,1) double {mustBeGreaterThanOrEqual(options.ChoiTolerance,0)}= 1e-12;
        options.decoySolver (1,1) string = "Mosek";
        options.decoyPrecision (1,1) string = "high";
        options.photonCutOff (1,1) double {mustBeInteger,mustBePositive} = 10;
        options.verboseLevel (1,1) double {mustBeInteger,mustBeNonnegative} = 1;
    end
   
    %number of intensities used
    n_decoy=numel(decoys);

    %Dimension of Bob
    dimB_SDP = size(POVMsBob{1},1);

    %Dimension of Alice
    dimA_SDP = size(statesAlice{1},1);

    %m = #number of rows of observations
    %n = #number of colums of observations
    m = size(expecJointTest,1);
    n = size(expecJointTest,2);
    
    %Empty vectors for upper and lower bounds
    Y1L = zeros(2,1);
    Y1U = ones(2,1);
    
    %cut-off for decoy
    n_photon= options.photonCutOff;
    
    % %Possonian distribution
    probDist=@(intensity,n) exp(-intensity).*intensity.^n./factorial(n);

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

    %Signal density matrices
    rhoAlice = cell(1,numel(statesAlice));
    for index = 1:numel(statesAlice)
        rhoAlice{index} = statesAlice{index}*statesAlice{index}';
    end

    %probability of sending in RL
    probRLCondDecoy = sum(expecJointTest(2:3,:,:),[1,2]);

    %QBER conditioned on each intensity
    qberAndDecoy = squeeze(expecJointTest(2,4,:) + expecJointTest(3,3,:));
    tqber = squeeze(t(2,4,:) + t(3,3,:));
    muqber = squeeze(mu(2,4,:) + mu(3,3,:));

    %Gain conditioned on each intensity
    gainAndDecoy = qberAndDecoy + squeeze(expecJointTest(2,3,:) + expecJointTest(3,4,:));
    tgain = tqber + squeeze(t(2,3,:) + t(3,4,:));
    mugain = muqber + squeeze(mu(2,3,:) + mu(3,4,:));

    %Calculate upper and lower bounds on the expectations
    upperbndsQber = min(max((qberAndDecoy + tqber + muqber)./(probsDecoy.*sum(probListAlice(2:3)).*ptest),0),1);
    lowerbndsQber = min(max((qberAndDecoy - tqber - muqber)./(probsDecoy.*sum(probListAlice(2:3)).*ptest),0),1);
    
    upperbndsGain = min(max((gainAndDecoy + tgain + mugain)./(probsDecoy.*sum(probListAlice(2:3)).*ptest),0),1);
    lowerbndsGain = min(max((gainAndDecoy - tgain - mugain)./(probsDecoy.*sum(probListAlice(2:3)).*ptest),0),1);
    
    %Projection onto n<=N_B subspace
    Proj_sub = blkdiag(eye(3),zeros(dimB_SDP-3));
    
    %solve for upper bound 1-photon    
    for indexopt = 1:1:2
            try
                %Set up LP for upper bounds
                % Y is the matrix representing all yields
                % J0 is the Choi matrix for the 0-photon component
                % J1 is the Choi matrix for the 1-photon component
                % In both Choi matrices Bob's system is first and Alice's
                % system is second to satisfy e.g. Y_0 = Tr[J_0 (F_y otimes (rho_x)^T)] 
                cvx_begin quiet
                    %CVX solver
                    cvx_solver(convertStringsToChars(options.decoySolver));
                    variable Y(n_photon+1,1)
                    variable J0(dimB_SDP,dimB_SDP) hermitian semidefinite
                    variable J1(dimA_SDP*dimB_SDP,dimA_SDP*dimB_SDP) hermitian semidefinite
                    maximize Y(2) 
                    subject to
                        % 0 <= Y <=1
                        vec(Y) >= vec(zeros(n_photon+1,1));
                        vec(Y) <= vec(ones(n_photon+1,1));
                        
                        %0- and 1-photon component treated seperately with add.
                        %constraints
                        Y0 = Y(1);
                        Y1 = Y(2);
                        
                        %Partial trace of Choi matrix = id
                        norm(PartialTrace(J1,[1],[dimB_SDP,dimA_SDP])-eye(dimA_SDP)) <= choi_tolerance;
                        norm(trace(J0) - 1) <= choi_tolerance;
                   
                        %Additional constraints for 0- and 1-photon components
                        %in terms of Choi matrices
                        if indexopt == 1 %QBER
                            %Usual decoy bounds rewritten as matrix times vector
                            for k = 1:n_decoy
                                pmu = probDist(decoys(k),0:1:n_photon);
                                vec(pmu*Y) >= vec(lowerbndsQber(k) - (1-sum(pmu)) - decoy_tolerance);
                                vec(pmu*Y) <= vec(upperbndsQber(k) + decoy_tolerance);                            
                            end
                            Y0 <= 1/2;
                            qberchannel1 = trace(J1*kron(POVMsBob{4},transpose(rhoAlice{2}))) + trace(J1*kron(POVMsBob{3},transpose(rhoAlice{3})));
                            norm(Y1 - qberchannel1) <= choi_tolerance;

                            qberchannel0 = trace(J0*kron(POVMsBob{4},1)) + trace(J0*kron(POVMsBob{3},1));
                            norm(Y0 - qberchannel0) <= choi_tolerance;
                        else %Gain
                            %Usual decoy bounds rewritten as matrix times vector
                            for k = 1:n_decoy
                                pmu = probDist(decoys(k),0:1:n_photon);
                                vec(pmu*Y) >= vec(lowerbndsGain(k) - (1-sum(pmu)) - decoy_tolerance);
                                vec(pmu*Y) <= vec(upperbndsGain(k) + decoy_tolerance);                            
                            end
                            gainchannel1 = trace(J1*kron(POVMsBob{3},transpose(rhoAlice{2}))) + trace(J1*kron(POVMsBob{4},transpose(rhoAlice{2}))) ...
                                            + trace(J1*kron(POVMsBob{3},transpose(rhoAlice{3}))) + trace(J1*kron(POVMsBob{4},transpose(rhoAlice{3})));
                            norm(Y1 - gainchannel1) <= choi_tolerance;

                            gainchannel0 = 2*(trace(J0*kron(POVMsBob{3},1)) + trace(J0*kron(POVMsBob{4},1)));
                            norm(Y0 - gainchannel0) <= choi_tolerance;
                        end
    
                        %Upper and lower bounds for n<=N_B subspace
                        for indexrow = 1:m
                            %1-photon
                            trace(J1*kron(Proj_sub,transpose(rhoAlice{indexrow}))) >= photonBound(indexrow,2) - choi_tolerance;                       
                            trace(J1*kron(Proj_sub,transpose(rhoAlice{indexrow}))) <= 1 + choi_tolerance;
                            %0-photon
                            trace(J0*kron(Proj_sub,1)) >= photonBound(indexrow,1) - choi_tolerance;
                            trace(J0*kron(Proj_sub,1)) <= 1 + choi_tolerance;
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
            Y1U(indexopt) = min(max(Y(2),0),1);            
    end
        
    % %solve for lower bounds
    for indexopt = 1:1:2
            try
                %Set up LP for upper bounds
                % Y is the matrix representing all yields
                % J0 is the Choi matrix for the 0-photon component
                % J1 is the Choi matrix for the 1-photon component
                % In both Choi matrices Bob's system is first and Alice's
                % system is second to satisfy e.g. Y_0 = Tr[J_0 (F_y otimes (rho_x)^T)] 
                cvx_begin quiet
                    %CVX solver
                    cvx_solver(convertStringsToChars(options.decoySolver));
                    variable Y(n_photon+1,1)
                    variable J0(dimB_SDP,dimB_SDP) hermitian semidefinite
                    variable J1(dimA_SDP*dimB_SDP,dimA_SDP*dimB_SDP) hermitian semidefinite
                    minimize Y(2) 
                    subject to
                        % 0 <= Y <=1
                        vec(Y) >= vec(zeros(n_photon+1,1));
                        vec(Y) <= vec(ones(n_photon+1,1));
                        
                        %0- and 1-photon component treated seperately with add.
                        %constraints
                        Y0 = Y(1);
                        Y1 = Y(2);
                        
                        %Partial trace of Choi matrix = id
                        norm(PartialTrace(J1,[1],[dimB_SDP,dimA_SDP])-eye(dimA_SDP)) <= choi_tolerance;
                        norm(trace(J0) - 1) <= choi_tolerance;
                   
                        %Additional constraints for 0- and 1-photon components
                        %in terms of Choi matrices
                        if indexopt == 1 %QBER
                            %Usual decoy bounds rewritten as matrix times vector
                            for k = 1:n_decoy
                                pmu = probDist(decoys(k),0:1:n_photon);
                                vec(pmu*Y) >= vec(lowerbndsQber(k) - (1-sum(pmu)) - decoy_tolerance);
                                vec(pmu*Y) <= vec(upperbndsQber(k) + decoy_tolerance);                            
                            end
                            Y0 <= 1/2;
                            qberchannel1 = trace(J1*kron(POVMsBob{4},transpose(rhoAlice{2}))) + trace(J1*kron(POVMsBob{3},transpose(rhoAlice{3})));
                            norm(Y1 - qberchannel1) <= choi_tolerance;
    
                            qberchannel0 = trace(J0*kron(POVMsBob{4},1)) + trace(J0*kron(POVMsBob{3},1));
                            norm(Y0 - qberchannel0) <= choi_tolerance;
                        else %Gain
                            %Usual decoy bounds rewritten as matrix times vector
                            for k = 1:n_decoy
                                pmu = probDist(decoys(k),0:1:n_photon);
                                vec(pmu*Y) >= vec(lowerbndsGain(k) - (1-sum(pmu)) - decoy_tolerance);
                                vec(pmu*Y) <= vec(upperbndsGain(k) + decoy_tolerance);                            
                            end
                            gainchannel1 = trace(J1*kron(POVMsBob{3},transpose(rhoAlice{2}))) + trace(J1*kron(POVMsBob{4},transpose(rhoAlice{2}))) ...
                                            + trace(J1*kron(POVMsBob{3},transpose(rhoAlice{3}))) + trace(J1*kron(POVMsBob{4},transpose(rhoAlice{3})));
                            norm(Y1 - gainchannel1) <= choi_tolerance;
    
                            gainchannel0 = 2*(trace(J0*kron(POVMsBob{3},1)) + trace(J0*kron(POVMsBob{4},1)));
                            norm(Y0 - gainchannel0) <= choi_tolerance;
                        end
    
                        %Upper and lower bounds for n<=N_B subspace
                        for indexrow = 1:m
                            %1-photon
                            trace(J1*kron(Proj_sub,transpose(rhoAlice{indexrow}))) >= photonBound(indexrow,2) - choi_tolerance;                       
                            trace(J1*kron(Proj_sub,transpose(rhoAlice{indexrow}))) <= 1 + choi_tolerance;
                            %0-photon
                            trace(J0*kron(Proj_sub,1)) >= photonBound(indexrow,1) - choi_tolerance;
                            trace(J0*kron(Proj_sub,1)) <= 1 + choi_tolerance;
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
            Y1L(indexopt) = min(max(Y(2),0),1);            
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