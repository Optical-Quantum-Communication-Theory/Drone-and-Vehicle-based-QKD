function [keyRate, modParser] = ThreeStateDataCarDifInt_KeyRateFunc(params,options,mathSolverFunc,mathSolverOptions,debugInfo)
% WLOG, this function works with the assumption that the z basis is used
% exclusively for generation rounds and the x basis is used exclusively for
% parameter estimation in the finite domain. Currently hard coded for two
% basis choices z and x.
% Aodhan Corrigan
%% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * announcementsA: Array of announcements made for each measurement Alice
%   made ordered in the same way as the columns of expectationsJoint.
% * announcementsB: Array of announcements made for each measurement Bob
%   made ordered in the same way as the rows of expectationsJoint.
% * keyMap: An array of KeyMapElement objects that contain pairs of accepted
%   announcements and an array dictating the mapping of Alice's measurement
%   outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   postive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that perform the pinching map 
%   key on  G(\rho). These projection operators should sum to identity.
% * f: error correction effiency. Set to 1 means for Shannon limit. 
% * observablesJoint: The joint observables of Alice and Bob's
%   measurments. The observables must be hermitian and each must be the size 
%   dimA*dimB by dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity. 
% * expectationsJoint: The joint expectations (as an array) from Alice and
%   Bob's measurements that line up with it's corresponding observable in
%   observablesJoint. These values should be betwen 0 and 1.
% * rhoA (nan): The fixed known density matrix on Alice's side for
%   prepare-and-measure protocols.
%% Following parameters are specific to Finite protocol %%
% * epsBounds: 1x4 vector containing epsilon bounds on subprotocols. Sum of
% all four elements should sum to total epsilon security. See George et al
% 2021 for details. 
% * t: Variation threshold on difference between observed statistical
% distribution and ideal statistical distribution \overline{F}. 
% * n, m: Signals used for key generation and parameter estimation
% respectively. Total signals sent is N = n + m.
% * Sigma, Lambda: Alphabet sizes. If no coarse graining is done, Sigma =
% Lambda. 
%% Outputs:
% * keyrate: Key rate of the QKD protocol measured in bits per block
%   processed.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * 
% 
% Reviewed by Devashish Tupkary 2023/09/18
% See also QKDKeyRateModule, FiniteBB84_4DAliceDescriptionFunc, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    mathSolverOptions (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;


%% modParser
modParser = moduleParser(mfilename);

%modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsJoint"]);

modParser.addRequiredParam("krausOps", @isCPTNIKrausOps);
modParser.addRequiredParam("keyProj", @(x) mustBeAKeyProj(x));
modParser.addRequiredParam("dimR",@mustBeInteger);

modParser.addRequiredParam("dimA",@mustBeInteger);
modParser.addRequiredParam("dimB", @mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"dimA")
modParser.addAdditionalConstraint(@mustBePositive,"dimB")
modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));

modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJoint","dimA","dimB"])

modParser.addRequiredParam("announcementsA")
modParser.addRequiredParam("genAnnouncementsA")
modParser.addRequiredParam("announcementsB")
modParser.addRequiredParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))
%modParser.addAdditionalConstraint(@mustBeSizedLikeAnnouncements,["expectationsJoint","announcementsA","announcementsB"])

%Epsilon(s)
modParser.addRequiredParam("epsilons", @(x) length(x) == 4);

%Error correction
modParser.addRequiredParam("f", @(x) mustBeGreaterThanOrEqual(x,1));

% Finite-size
% modParser.addRequiredParam("n", @(n) mustBeInteger(n));
% modParser.addRequiredParam("m", @(m) mustBeInteger(m));
modParser.addRequiredParam("N", @(N) mustBeInteger(N));
modParser.addRequiredParam("t", @mustBeNonnegative);
modParser.addRequiredParam("tsift", @mustBeNonnegative);

% modParser.addRequiredParam("pzB", @(pz) mustBeInRange(pz, 0, 1));
modParser.addRequiredParam("pGen", @(p) mustBeInRange(p, 0, 1)); 
modParser.addOptionalParam("rhoA", nan, @(x) all(isnan(x),"all") ||isDensityOperator(x));
modParser.addOptionalParam("rhoAAp", nan, @(x) all(isnan(x),"all") ||isDensityOperator(x));

modParser.addOptionalParam("blockDimsA", nan); % figure out check on these later (both are given or not, dimensions sum to total dimensions etc.)
modParser.addOptionalParam("blockDimsB", nan);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsA","dimA"]);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsB","dimB"]);
modParser.addAdditionalConstraint(@(blockDimsA,blockDimsB) ~xor(all(isnan(blockDimsA),"all"),all(isnan(blockDimsB),"all")),["blockDimsA","blockDimsB"]);

% modParser.addRequiredParam("Sigma", @(x) mustBeInteger(x));
% modParser.addRequiredParam("Lambda", @(x) mustBeInteger(x));
% modParser.addOptionalParam("coarseGrain", @(x) size(x,1) == params.Sigma || size(x,2) == params.Lambda);

% modParser.addRequiredParam("expectationsConditional",@(x) eachRowMustBeAProbDist(x)); %%P( b | a, mu_i, test *or* gen, coherent state sent into channel)
% modParser.addRequiredParam("decoys",@(x) allCells(x,@(y) y>=0));

%Signal intensity
modParser.addRequiredParam("decoysSignal",@(x) allCells(x,@(y) y>=0));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoysSignal");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoysSignal");
%Decoy intensity 1
modParser.addRequiredParam("decoys1",@(x) allCells(x,@(y) y>=0));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoys1");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoys1");
%Decoy intensity 2
modParser.addRequiredParam("decoys2",@(x) allCells(x,@(y) y>=0));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoys2");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoys2");

modParser.addRequiredParam("decoyProbs",@mustBeProbDistCell); %,must sum to 1. 
modParser.addRequiredParam("probSignalsAgen",@mustBeProbDist); %Add something for same number of announcementsA or something
modParser.addRequiredParam("probSignalsAtest",@mustBeProbDist); %Add something for same number of announcementsA or something
modParser.addRequiredParam("signalsAlice", @(x) length(x) == 3); 
modParser.addRequiredParam("POVMB", @(x) length(x) == 8); 
%modParser.addRequiredParam("observablesJointTesting");
%modParser.addRequiredParam("observablesJointGeneration");
modParser.addRequiredParam("pTestandACondSingle");
modParser.addRequiredParam("pMuandACondTest");
modParser.addRequiredParam("pSinglePhoton", @(p) p >= 0);
modParser.addRequiredParam("pGenCondSingle", @(p) p >= 0); 
modParser.addRequiredParam("pACondSingle", @mustBeProbDist); 
modParser.addRequiredParam("patternOrder"); 
modParser.addRequiredParam("pTestCondSingleandA");
modParser.addRequiredParam("pSingleCondGen");

modParser.addRequiredParam("tau", @(tau) all(size(tau) == [3, 8, 3])); %Alice & Bob alphabets, mu intensities

modParser.addRequiredParam("ohioData"); 
modParser.addRequiredParam("ohioOrder"); 
modParser.addRequiredParam("ohioObservations"); 

modParser.addRequiredParam("cn"); 


modParser.parse(params);
params = modParser.Results;

optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.addOptionalParam("decoyTolerance",1e-14,@(x) x>=0);
optionsParser.addOptionalParam("decoySolver","Mosek");
optionsParser.addOptionalParam("decoyPrecision","default");
optionsParser.parse(options);
options = optionsParser.Results;

%% simple setup
debugMathSolver = debugInfo.addLeaves("mathSolver");
mathSolverInput = struct();

%% Import Data
pGen = params.pGen; 
decoyProbs = cell2mat(params.decoyProbs); 
probSignalsTesting = params.probSignalsAtest;
probSignalsGeneration = params.probSignalsAgen; 
pMuandACondTest = params.pMuandACondTest;

% Define squashing map
squashingMap = ThreeStateFlagStateSquashingPostProccessingMap();

% Convert data into our format
data = convertData(params.ohioObservations, params.ohioOrder, params.patternOrder);

%apply squashing map

%P(b | a, gen) 
squashedCondExpSignalOnly = data(1:3, :)*squashingMap.'; 

%P(b| a , mu_i, test) 
squashedConExp(:, :, 1) = data(4:6, :)*squashingMap.';
squashedConExp(:, :, 2) = data(7:9, :)*squashingMap.';
squashedConExp(:, :, 3) = data(10:12, :)*squashingMap.';

%start P(b | a, mu, test) NO PHOTON CONDITIONING
squashedExpCondionedOnTest = zeros(size(squashedConExp)); %initialize for storage

%calculate unconditional expectations
for i=1:numel(decoyProbs)
    squashedExpCondionedOnTest(:,:,i) = diag(pMuandACondTest(:, i))*squashedConExp(:,:,i);% P(a,b,mu | test)
end

%P(a,b,mu,test)
squashedJointExp = (1-pGen)*squashedExpCondionedOnTest; %P(a,b,mu,test)



%% Calculate mu ball
obsTau = params.tau; 

epsAT = params.epsilons(4);

[muLower,muUpper] = mubetaUpperLower(params.N, squashedJointExp, params.tau, log10(epsAT));
muMat = max(muUpper,muLower);

muBall = min(muMat,[],"all");
params.muBall = muBall;
debugInfo.storeInfo("mu", muBall);
%muBall = 0;

%% Error correction
squashedExpCondMuGen = diag(probSignalsGeneration)*squashedCondExpSignalOnly; %P(a,b|gen)
[deltaLeak, gains] = errorCorrectionCost(params.genAnnouncementsA,params.announcementsB,...
    squashedExpCondMuGen,params.keyMap,params.f); 

debugInfo.storeInfo("deltaLeak",deltaLeak);
debugInfo.storeInfo("gains", gains);

%% translate for the math solver
%now we have all the parts we need to get a key rate from the a math
%solver, but we need to put it into a form it can understand.
%first we give it the kraus operators for the G map and the projection
%operators for the key map (Z).
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;
% also include rhoA from the description if it was given

if ~isnan(params.rhoA)
    mathSolverInput.rhoA = params.rhoA;
end

%% Decoy Analysis %%
reshapedIntensities = {cell2mat(params.decoysSignal);cell2mat(params.decoys1);cell2mat(params.decoys2)};
decoyIntensityFlag = reshapedIntensities{2};
pTestCondSingleandA = params.pTestCondSingleandA; 

%Generate photon bounds for decoy analysis
%observations for flag-state squasher used in Decoy SDP
obsFlagDecoy = 1 - sum(squashedConExp(:,1:4,:) - obsTau(:,1:4,:) - muMat(:,1:4,:),2) - squashedConExp(:,end,:) + obsTau(:,end,:) + muMat(:,end,:);
for indexA = 1:numel(params.signalsAlice)
    photonBoundDecoy(:,2) = 1 - obsFlagDecoy(:,:,2)./(decoyIntensityFlag(indexA)*exp(-decoyIntensityFlag(indexA))*params.cn); % 1-photon
    photonBoundDecoy(:,1) = 1 - obsFlagDecoy(:,:,2)./(exp(-decoyIntensityFlag(indexA))*params.cn); % 0-photon
end

%Fine-grained analysis

% % SDP
% [CondExpectationsL, CondExpectationsU] = decoySDPFiniteDifInt(squashedJointExp,params.signalsAlice,params.POVMB,reshapedIntensities',probSignalsTesting,decoyProbs,1-pGen,photonBoundDecoy,params.tau,muMat, ...
%    "decoySolver",options.decoySolver,"decoyTolerance",options.decoyTolerance,...
%     "decoyPrecision",options.decoyPrecision,"verboseLevel",options.verboseLevel);

%LP
[CondExpectationsL, CondExpectationsU] = decoyLPFiniteDifInt(squashedJointExp,reshapedIntensities',probSignalsTesting,decoyProbs,1-pGen,params.tau,muMat, ...
   "decoySolver",options.decoySolver,"decoyTolerance",options.decoyTolerance,...
    "decoyPrecision",options.decoyPrecision,"verboseLevel",options.verboseLevel);

% %Independent LP
% [CondExpectationsL, CondExpectationsU] = decoyIndependent(squashedJointExp,reshapedIntensities,probSignalsTesting,decoyProbs,1-pGen,params.tau,muMat, ...
%    "decoySolver",options.decoySolver,"decoyTolerance",options.decoyTolerance,...
%     "decoyPrecision",options.decoyPrecision,"verboseLevel",options.verboseLevel);


% We convert them to Prob (Bob and Alice | one-photon, key/test, intensity);
JointExpectationsL = diag(params.pTestandACondSingle)*CondExpectationsL; %P(a,t|n)P(b|a,t,n)=P(a,b,t|n)
JointExpectationsU = diag(params.pTestandACondSingle)*CondExpectationsU;

%Now, we add constraints. Recall that the POVMs correspond to Bob and Alice
%given Test . (and trivially intensity, and single photon)

numObs = numel(params.observablesJoint);
observablesJointTesting = params.observablesJoint; 
for a = 1:3
    for b = 1:8
        observablesJointTesting{a,b} = pTestCondSingleandA(a) * observablesJointTesting{a,b};
    end
end

% Decoy Inequality Constraints
decoyConstraints = arrayfun(@(index)InequalityConstraint(...
    observablesJointTesting{index},JointExpectationsL(index),...
    JointExpectationsU(index)), 1:numObs);

%% Generate photon number constraints for a = 1, 2, 3 from decoy analysis
pU = params.pACondSingle; 
pL = zeros([3,1]);

obsFlagKeyOpt = 1 - sum(squashedConExp(:,1:4,:) - obsTau(:,1:4,:) - muMat(:,1:4,:),2) - squashedConExp(:,end,:) + obsTau(:,end,:) + muMat(:,end,:);
decoyIntensityFlagKeyOpt = reshapedIntensities{1};
for indexA =1:numel(params.signalsAlice)
    photonBoundKeyOpt = 1 - obsFlagKeyOpt(:,:,1)./(decoyIntensityFlagKeyOpt(indexA)*exp(-decoyIntensityFlagKeyOpt(indexA))*params.cn); % 1-photon
end

for i = 1:3 
    % %CondExpectationsL  S_( H V D A ) D_( HV RL ) CC_(ANY) NON
    singleClicksM = sum(CondExpectationsL(i, 1:4)); 
    vacM = CondExpectationsL(i, end); 
    multiClicksM = sum(CondExpectationsU(i, 5:end-1)); 
    mObs = min(1 - singleClicksM - vacM, multiClicksM);
    %fprintf('%f\n', params.cn);
    pL(i) = params.pACondSingle(i)* (1 - mObs/(params.cn)); 
    % pL(i) = params.pACondSingle(i)*photonBoundKeyOpt(i); 
end

%Constraint ops to project onto n <= NB subspace
projectionOps = cell([1, 3]);
Pi = blkdiag(eye(3), zeros(7)); 
for i = 1:3
    aa = zket(3, i)*zket(3, i)';
    projectionOps{i} = kron(aa, Pi); 
end
photonNumberConstraints = arrayfun(@(index) InequalityConstraint(projectionOps{index}, pL(index), pU(index)), 1:3); 


%% Define final inequality constraints
ineqConstraintsTotal = [decoyConstraints, photonNumberConstraints]; 

%Add inequality constraints to math solver
mathSolverInput.inequalityConstraints = ineqConstraintsTotal; 


% if block diag information was give, then pass it to the solver.
if ~any(isnan(params.blockDimsA),"all")
    mathSolverInput.blockDimsA = params.blockDimsA;
    mathSolverInput.blockDimsB = params.blockDimsB;
end

% now we call the math solver function on the formulated inputs, with
% it's options.
[relEnt,~] = mathSolverFunc(mathSolverInput,mathSolverOptions,debugMathSolver);
debugInfo.storeInfo("relEnt",relEnt);

%store the key rate (even if negative)
% [keyRate,ECCost,AEPCorrectionPD,privacyAmp] = computeFiniteKeyRate(relEnt, deltaLeak, gains, obsTau, params);
[keyRate,ECCost,AEPCorrectionPD,privacyAmp,debugInfo] = finiteKeyRate(relEnt, deltaLeak, gains, params, options, debugInfo);

debugInfo.storeInfo("ECCost", ECCost); 
debugInfo.storeInfo("AEPCorrectionPD", AEPCorrectionPD); 
debugInfo.storeInfo("privacyAmp", privacyAmp); 

if options.verboseLevel>=1
    %ensure that we cut off at 0 when we display this for the user.
    fprintf("Key rate: %e\n",max(keyRate,0));
end

%set the linearization estimate key rate as well for debuging
if isfield(debugMathSolver.info,"relEntStep2Linearization")
    keyRateStep2Linearization = debugMathSolver.info.relEntStep2Linearization - deltaLeak ; 
    debugInfo.storeInfo("keyRaterelEntStep2Linearization",keyRateStep2Linearization)

    if options.verboseLevel>=2
        fprintf("Key rate using step 2 linearization intial value: %e\n",max(keyRateStep2Linearization,0))
    end
end
end

function [keyRate, ECLeakage, AEPCorrectionPD, privacyAmplification, debugInfo] = finiteKeyRate(relEnt, deltaLeak, gains, params, options, debugInfo)
    %computes the finite size keyrate.
    
    %Alice's key alphabet size
    dA = params.dimR;
    
    %key generation probability
    pgen = params.pGen;
    
    %total signals sent
    N = params.N;

    %epsilons
    epsBar = params.epsilons(1); 
    epsEC = params.epsilons(2);
    epsPA = params.epsilons(3);
    epsAT = params.epsilons(4);
    
    %t parameter for allowed fluctuations
    t = params.t;
    tsift = params.tsift;
    
    %AEP correction in purified distance
    AEPCorrectionPD = 2*log2(2*dA+1)*sqrt(1-2*log2(epsBar));
    
    %Privacy Amplification
    privacyAmplification = 2*log2(1/2/epsPA)/N;
    
    %Error correction including error verification
    % deltaleak calculated conditioned on having a generation round, hence the
    % prefactor of pgen
    ECLeakage = pgen*deltaLeak + log2(2/epsEC)/N; %actually epsilon.EV this is!
    
    %Total gain
    totalGain = sum(gains,1);
    
    %Sifted number of signals
    n_sift = pgen*totalGain*N; 
    
    %Mu sift
    [~,musift] = mubetaUpperLower(N,totalGain*pgen,tsift,log10(epsAT));
    
    %Probabilty of Alice sending a single photon 
    % pSinglePhoton = params.decoys{1}*exp(-params.decoys{1});
    
    %probability of a single photon in a generation round
    pSinglePhotonCondGen = params.pSingleCondGen;
    pGenCondSingle = params.pGenCondSingle;
    
    %Calculate prefactor of relative entropy
    % prefactor = pgen*pSinglePhoton/N*floor(n_sift-tsift*N)/(totalGain*pgen + tsift + musift); 
    prefactor = pSinglePhotonCondGen*floor(n_sift-tsift*N)/(sum(gains)*pGenCondSingle + tsift + musift)/N;
    
    %Final resulting finite-size key rate
    keyRate = prefactor*relEnt - sqrt((totalGain*pgen - tsift)/N)*AEPCorrectionPD - ECLeakage - privacyAmplification; 

end


function mapping = ThreeStateFlagStateSquashingPostProccessingMap()
% squashes detector patterns from H,V,R,L (as bit string patterns 0000,
% 1000, ..., 1111) to signal+vac values H,V,R,L, (double clicks), CC, vac.

function mapping = quickMap(mapping,pattern,remaping)
    mapping(:,sub2indPlus(2*ones(1,numel(pattern)),pattern+1)) = remaping;
end

mapping = zeros(8,16);

% A^T = Lambda * E^T => A = E*Lambda^T
% => Lambda 8 by 16 since E 3x16 and A 3x8

% The vast majority of the squashed bits are cross clicks that are mapped
% to vac for discarding. We will replace patterns that don't represent
% cross clicks in later steps.
%mapping(5,:) = 1;

%order S_( H V D A ) D_( HV RL ) CC_(ANY) NON

% zero clicks to vac
mapping = quickMap(mapping,[0,0,0,0],[0,0,0,0,0,0,0,1]);

% single clicks to single clicks
mapping = quickMap(mapping,[1,0,0,0],[1,0,0,0,0,0,0,0]); % H
mapping = quickMap(mapping,[0,1,0,0],[0,1,0,0,0,0,0,0]); % V
mapping = quickMap(mapping,[0,0,1,0],[0,0,1,0,0,0,0,0]); % R
mapping = quickMap(mapping,[0,0,0,1],[0,0,0,1,0,0,0,0]); % L

% double clicks
mapping = quickMap(mapping,[1,1,0,0],[0,0,0,0,1,0,0,0]); % Z (HV)
mapping = quickMap(mapping,[0,0,1,1],[0,0,0,0,0,1,0,0]); % X (RL)

%cross clicks
% 2 clicks
mapping = quickMap(mapping,[1,0,1,0],[0,0,0,0,0,0,1,0]); %HR
mapping = quickMap(mapping,[1,0,0,1],[0,0,0,0,0,0,1,0]); %HL
mapping = quickMap(mapping,[0,1,1,0],[0,0,0,0,0,0,1,0]); %VR
mapping = quickMap(mapping,[0,1,0,1],[0,0,0,0,0,0,1,0]); %VL
% 3 clicks
mapping = quickMap(mapping,[1,1,1,0],[0,0,0,0,0,0,1,0]); %HV R
mapping = quickMap(mapping,[1,1,0,1],[0,0,0,0,0,0,1,0]); %HV L
mapping = quickMap(mapping,[1,0,1,1],[0,0,0,0,0,0,1,0]); %H RL
mapping = quickMap(mapping,[0,1,1,1],[0,0,0,0,0,0,1,0]); %V RL
%4 clicks
mapping = quickMap(mapping,[1,1,1,1],[0,0,0,0,0,0,1,0]); %HV RL

end


%Checker that each row is a probability distribution
function eachRowMustBeAProbDist(expectationsConditional)

% get the dimensions of the conditional expectations. Then based on that
% pick a strategy to handle it
dimExpCon = size(expectationsConditional);

errorID ="ThreeStateDecoyKeyRateFunc:InvalidRowsAreNotProbDists";
errorTXT = "A row in the conditional distribution is not a valid probability distribution.";

if numel(dimExpCon) == 2 % Matlab's minimum number of dimensions is 2.
    % The array is 2d and the slicing is easy
    for index = 1:dimExpCon(1)
        if~isProbDist(expectationsConditional(index,:))
           throwAsCaller(MException(errorID,errorTXT));
        end
    end
else
    % We have some tricky slicing to do for 3 plus dimensions.
    % We need to index the first dimension and the combination of
    % dimensions 3 and up. The second dimension will just use :.
    maskedDims = [dimExpCon(1),prod(dimExpCon(3:end))];

    for index = 1:prod(maskedDims)
        vecIndex = ind2subPlus(maskedDims,index);
        if ~isProbDist(expectationsConditional(vecIndex(1),:,vecIndex(2)))
            throwAsCaller(MException(errorID,errorTXT));
        end
    end
end
end

function observablesAndDimensionsMustBeTheSame(observables,dimA,dimB)
if ~allCells(observables,@(x) size(x,1) == dimA*dimB)
    throwAsCaller(MException("FiniteKeyRateFunc:ObservablesAndDimensionsMustBeTheSame","The Observables must have the same dimensions as Alice and Bob multiplied together."));
end
end
function mustBeProbDistCell(input)
mustBeProbDist([input{:}])
end
function mustBeSizedLikeAnnouncements(jointExpectations,announcementsA,announcementsB)
if ~isequal(size(jointExpectations),[numel(announcementsA),numel(announcementsB)])
    throwAsCaller(MException("FiniteKeyRateFunc:jointKeyDoesNotHaveSameSizeAsAnnouncements",...
        "The joint key distribution must have size numel(announcementsA) by numel(announcementsB)."))
end
end
