%helper function that performs decoy state analysis
function [Y1L,Y1U] = decoyAnalysisConstrFiniteOhio(expecJointTest,decoys,probListAlice,probsDecoy,ptest,mu,options)
    arguments
        expecJointTest (:,:,:) double {mustBeGreaterThanOrEqual(expecJointTest,0)}
        decoys (:,1) double {mustBeGreaterThanOrEqual(decoys,0),SameNumberOfDecoysAndPages(decoys,expecJointTest)}        
        probListAlice (:,1) double {mustBeProbDist}
        probsDecoy (:,1) double {mustBeProbDist}
        ptest (1,1) double {mustBeInRange(ptest,0,1),mustBeSameTotalProb(ptest,expecJointTest)}
        mu (1,1) double {mustBeGreaterThanOrEqual(mu,0)}
        options.t (1,1) double = nan
        options.tau (:, :, :) double = nan
        options.decoyTolerance (1,1) double {mustBeGreaterThanOrEqual(options.decoyTolerance,0)}= 1e-14;
        options.decoySolver (1,1) string = "Mosek";
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
    Y1U = zeros(m,n);
    
    %cut-off for decoy
    n_photon= options.photonCutOff;
    
    %Possonian distribution
    Poisson=@(intensity,n) exp(-intensity).*intensity.^n./factorial(n);
    
    %Tolerance in Decoy LP
    decoy_tolerance = options.decoyTolerance;

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
    upperbndsExpec = zeros(m,n,numel(decoys));
    lowerbndsExpec = zeros(m,n,numel(decoys));
    
   if isnan(options.t) == isnan(options.tau)
        %check both defined or neither defined
        throw(MException('decoyAnalysisConstrFinite:tVals', ...
        'Must specify exactly one of t or tau.'));
   elseif ~isnan(options.tau)
        %if tau defined
        try
            mustBeGreaterThanOrEqual(options.tau, 0);
            mustBeCompatibleSize(options.tau, expecJointTest);
            %check >0 and compatible size
        catch error
            rethrow(error);
        end
        tau = options.tau; 
   elseif ~isnan(options.t)
       %if t defined
        try
            mustBeGreaterThanOrEqual(options.t, 0);
            %check >0
        catch error
            rethrow(error);
        end
        tau = options.t*ones(size(expecJointTest)); 
   end
    for kdecoy = 1:numel(decoys)
        for row = 1:numel(probListAlice)
            upperbndsExpec(row,:,kdecoy) = min(max((expecJointTest(row,:,kdecoy) + tau(row, :, kdecoy) + mu)/(probListAlice(row)*probsDecoy(kdecoy)*ptest),0),1);
            lowerbndsExpec(row,:,kdecoy) = min(max((expecJointTest(row,:,kdecoy) - tau(row, :, kdecoy) - mu)/(probListAlice(row)*probsDecoy(kdecoy)*ptest),0),1);
        end
    end
    
    %solve for upper bound 1-photon    
    for i=0:m-1
        for j=1:n

            if options.verboseLevel>=2
                fprintf("\n Upper bounds 1-photon solving LP %d of %d ",n*i+(j-1)+1,m*n)
            end
            
            try
                %Set up LP for upper bounds
                % Y is the matrix representing all yields
                cvx_begin quiet
                    %CVX solver
                    cvx_solver(convertStringsToChars(options.decoySolver));
                    variable Y(m*(n_photon+1),n)
                    maximize Y(i*(n_photon+1)+2,j)
                    subject to
                        % 0 <= Y <=1       
                        Y >= 0;
                        Y <= 1;
                        
                        %{
                        %0-photon yields
                        Y0 = M{1}*Y;
                        
                        %0-photon error rate e_0=1/2
                        Y0(1,2) + Y0(2,1) == 1/2*(Y0(1,1)+Y0(1,2)+Y0(2,1)+Y0(2,2));
                        Y0(3,4) + Y0(4,3) == 1/2*(Y0(3,3)+Y0(3,4)+Y0(4,3)+Y0(4,4));
                        %}
                        %Usual decoy bounds rewritten as matrix times vector
                        for k = 1:n_decoy
                            pmu = Poisson(decoys(k),0:1:n_photon);
                            Pmu = kron(eye(m),pmu);
                            Pmu*Y >= lowerbndsExpec(:,:,k) - ones(m,n) * (1-sum(pmu)) - decoy_tolerance;
                            Pmu*Y <= upperbndsExpec(:,:,k) + decoy_tolerance;                            
                        end
                cvx_end

                if options.verboseLevel>=2
                    disp(cvx_status)
                end

            if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
                Y(i*(n_photon+1)+2,j) = 1;
            end
            catch error 
                rethrow(error);%  fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end

            %Store upper bound
            Y1U(i+1,j) = Y(i*(n_photon+1)+2,j);            
        end
    end
        
    %solve for lower bound
    for i=0:m-1
        for j=1:n

            if options.verboseLevel>=2
                fprintf("\n Lower bounds 1-photon solving LP %d of %d ",n*i+(j-1)+1,m*n)
            end
            
            try
                %Set up LP for upper bounds
                % Y is the matrix representing all yields
                cvx_begin quiet
                    %CVX solver
                    cvx_solver(convertStringsToChars(options.decoySolver));
                    variable Y(m*(n_photon+1),n)
                    minimize Y(i*(n_photon+1)+2,j)
                    subject to

                        % 0 <= Y <=1
                        Y >= 0;
                        Y <= 1;
                        
                        %{
                        %0-photon yields
                        Y0 = M{1}*Y;

                        %0-photon error rate e_0=1/2
                        Y0(1,2) + Y0(2,1) == 1/2*(Y0(1,1)+Y0(1,2)+Y0(2,1)+Y0(2,2));
                        Y0(3,4) + Y0(4,3) == 1/2*(Y0(3,3)+Y0(3,4)+Y0(4,3)+Y0(4,4));
                        %}
                        %Usual decoy bounds rewritten as matrix times vector
                        for k = 1:n_decoy
                            pmu = Poisson(decoys(k),0:1:n_photon);
                            Pmu = kron(eye(m),pmu);
                            Pmu*Y >= lowerbndsExpec(:,:,k) - ones(m,n) * (1-sum(pmu)) - decoy_tolerance;
                            Pmu*Y <= upperbndsExpec(:,:,k) + decoy_tolerance;                            
                        end
                cvx_end

                if options.verboseLevel>=2
                    disp(cvx_status)
                end

            if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
                Y(i*(n_photon+1)+2,j) = 0;
            end
            catch 
                rethrow(error); %fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end

            %Store lower bound
            Y1L(i+1,j) = Y(i*(n_photon+1)+2,j);
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

function mustBeCompatibleSize(x,y)
    if  size(x) ~= size(y)
        throw(MException('mu_beta_func:FrequenciesTauMustMatch', ...
        'Dimension of the Frequency Array and t Values Array must match.'));
    end
end