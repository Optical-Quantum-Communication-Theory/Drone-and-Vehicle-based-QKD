function [newParams,modParser] = ThreeStateDataAirtoCarDifInt_Desc(params, options, debugInfo)
% BasicBB84_LossyDescriptionFunc A simple description function for a qubit BB84
% protocol with loss, using the Schmidt decomposition to turn Alice's
% 4d space of signals sent to a 2d space. 
%
% Input parameters:
% * pzA/B: The probability that Alice and Bob measure in the Z-basis. 
% Output parameters:
% * observablesJoint: The joint observables for Alice and Bob's measurement
%   of the signals.
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * probSignalsA: Probability of Alice selecting a signal to send to Bob.
%   In this protocol, it is half the probability of the basis choice.
% * rhoA: Alice's reduced density matrix for prepare-and-measure based
%   protocols.
% * POVMA: Alice's set of POVM operators which she measures her state with
%   in the source replacement scheme.
% * POVMB: Bob's set of POVM operators which he measures his state with.
% * announcementsA: Alice's announcements for each of her POVM operators.
%   Can be integers or strings.
% * announcementsB: Bob's announcements for each of his POVM operators.
%   Can be integers or strings.
% * keyMap: An array of KeyMap objects that contain pairs of accepted
%   announcements and an array dictating the mapping of Alice's measurement
%   outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   postive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that extract the key
%   from G(\rho). These projection operators should sum to identity. This
%   map is often called Z.
% Options:
% * none
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i which should be <= identity for
%   a CPTNI map.
%
% Reviewed by Devashish Tupkary 2023/09/18
% See also QKDDescriptionModule, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
%Parsing technical options for the module
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
%Parsing parameters for the module
modParser = moduleParser(mfilename);

%Finite size
modParser.addRequiredParam("N", @(N) mustBeInteger(N));
modParser.addRequiredParam("epsSound", @(e) mustBeNonnegative(e));
modParser.addRequiredParam("epsATfrac", @(f) mustBeInRange(f, 0, 1));
%modParser.addRequiredParam("t", @(t) mustBeNonnegative(t)); 
modParser.addRequiredParam("tau", @(t) mustBeNonnegative(t)); 

%Protocol parameters
modParser.addRequiredParam("pGen",@(x) mustBeInRange(x,0,1));
modParser.addRequiredParam("pzB",@(x) mustBeInRange(x,0,1));
modParser.addRequiredParam("decoysSignal",@(x) allCells(x,@(y) y>=0));
modParser.addRequiredParam("decoys1",@(x) allCells(x,@(y) y>=0));
modParser.addRequiredParam("decoys2",@(x) allCells(x,@(y) y>=0));
modParser.addRequiredParam("decoyProbs",@mustBeProbDistCell); %,must sum to 1. 

modParser.addRequiredParam("ohioG");

modParser.parse(params)
params = modParser.Results;

%% simple setup
%Setup parameters for later use
newParams = struct();

pGen = params.pGen; 
pTest = 1 - pGen; 

%Define R and L states for later use
ketR = [0.695053523375061; -0.124886899941761 + 0.708028150476271*1i];
ketL = [0.681617859201278; 0.140536423004987 - 0.718085376419009*1i];
ketH = [0.997915608957927; -0.015533489646278 - 0.062635038896280*1i];

%store Alice's sent states
newParams.signalsAlice = {ketH,ketR,ketL};

dimA = 3; % |1,2,3>
dimB = 10; % comb basis + 2x DD + CC + VAC flags

newParams.dimA = dimA;
newParams.dimB = dimB;

%Dimension of Alice's key register
newParams.dimR = 2;

pzB = params.pzB; 
pxB = 1 - pzB; 

%% generate rhoA
%% PROBABILITIES FROM EXPERIMENT %%
pAandGen = [0, (51/256)*(119/121), (51/256)*(119/121)]; 
pAandTest = [1107/2560, 54433/619520, 54433/619520];  
pGen = sum(pAandGen); 
pTest = sum(pAandTest); 
probSignalsAgen = pAandGen/pGen; %p(a|gen)
probSignalsAtest = pAandTest/pTest; %p(a|test)


newParams.probSignalsAgen = probSignalsAgen; 
newParams.probSignalsAtest = probSignalsAtest; 
newParams.pGen = pGen; 
newParams.pTest = pTest; 

%Store intensities in a matrix
reshapedIntensities = [cell2mat(params.decoysSignal);cell2mat(params.decoys1);cell2mat(params.decoys2)];


[pA, pTestCondSingle, pGenCondSingle, pTestCondSingleandA, pGenCondSingleandA, pACondTestandSingle, pACondSingle, pSinglePhoton ,pSingleCondGen,pMuandACondTest] = ...
    computeProbabilitiesFuncDifInt(pGen, pTest, probSignalsAgen, probSignalsAtest, reshapedIntensities, params.decoyProbs);

newParams.pSinglePhoton = pSinglePhoton; 
debugInfo.storeInfo("pSinglePhoton", pSinglePhoton); 
newParams.pTestandACondSingle = pTestCondSingle*pACondTestandSingle; %p(a,test|n=1)
newParams.pGenCondSingle = pGenCondSingle; %p(gen|n=1)
newParams.pACondSingle = pACondSingle; %p(a|n=1)
newParams.pTestCondSingleandA = pTestCondSingleandA; %p(test|n=1,a) for a =1,2,3
newParams.pMuandACondTest = pMuandACondTest; %p(mu,a|test)
debugInfo.storeInfo("pMuandACondTest", pMuandACondTest); 
newParams.pSingleCondGen = pSingleCondGen;

%condition on generation for setup

psi1 = sqrt(pACondSingle(1))*kron(zket(3,1), ketH);
psi2 = sqrt(pACondSingle(2))*kron(zket(3,2), ketR);
psi3 = sqrt(pACondSingle(3))*kron(zket(3,3), ketL);

psiAAp = psi1 + psi2 + psi3;  %% |psi> conditioned on n = 1. 
rhoAAp = psiAAp * psiAAp'; 
newParams.rhoAAp = rhoAAp; 

pert = perturbationChannelEpsilon(rhoAAp);
rhoAApPert = perturbationChannel(rhoAAp,pert);

rhoA = PartialTrace(rhoAAp, 2, [3, 2]); 
newParams.rhoA = rhoA; 

%% joint obserables
%testPOVMsA = {pTestCondSingleA(1)*diag([1,0,0]), pTestCondSingleA(2)*diag([0,1,0]), pTestCondSingleA(3)*diag([0,0,1])};
%genPOVMsA = {pGenCondSingleA(2)*diag([0,1,0]), pTestGenSingleA(3)*diag([0,0,1])};
POVMsA = {diag([1,0,0]), diag([0,1,0]), diag([0,0,1])}; 


%Flag state POVMsB
%order S_( H V R L ) D_( HV RL ) CC_(ANY) NON
G = params.ohioG; 
G = conj(G([3, 4, 1, 2, 7, 8, 5, 6], :)); 

[Vd, newParams.cn,~,~] = detectorDecompFlag(G, 1, false); 
leqOneSubpsacePOVMs = newBasisPOVMsB(Vd); 

GammaH = blkdiag(leqOneSubpsacePOVMs{1}, diag(zket(7,1)));
GammaV = blkdiag(leqOneSubpsacePOVMs{2}, diag(zket(7,2)));
GammaR = blkdiag(leqOneSubpsacePOVMs{3}, diag(zket(7,3)));
GammaL = blkdiag(leqOneSubpsacePOVMs{4}, diag(zket(7,4)));
GammaHV = blkdiag(zeros(3), diag(zket(7,5)));
GammaRL = blkdiag(zeros(3),diag(zket(7,6)));
GammaCC = blkdiag(zeros(3), diag(zket(7,7)));
GammaVAC = blkdiag(diag([0,0,1]), zeros(7));

POVMsB = {GammaH, GammaV, GammaR, GammaL, GammaHV, GammaRL, GammaCC, GammaVAC};
debugInfo.storeInfo('POVMsB', POVMsB);
mustBePOVM(POVMsB); 

newParams.POVMA = POVMsA;
newParams.POVMB = POVMsB;


% each POVM element is assigned an announcement made by their respective
% party
newParams.announcementsA = ["GZ", "TZ", "GX", "TX", "GX", "TX"];

newParams.genAnnouncementsA = ["GZ", "GX", "GX"];
newParams.announcementsB = ["Z","Z","X","X", "DZ", "DX", "CC", "VAC"];

newParams.keyMap = [ KeyMapElement("GX", "X", ["T", 1, 2]) ]; %% Z basis used for testing 
%First 1 should be N, assigned to 1 to save on memory since will never be
%used. 

observablesJoint = cell(numel(POVMsA), numel(POVMsB)); 
for a = 1:numel(POVMsA)
    for b = 1:numel(POVMsB)
        observablesJoint{a,b} = kron(POVMsA{a}, POVMsB{b}); 
    end
end     

newParams.observablesJoint = observablesJoint; 

%% Kraus Ops (for G map)

%K_X shrunk via iso = |0X2|_RA + |1X3|_RA. 
[P, D] = eig(GammaR+GammaL);
sqrtEigenvals = real(sqrt(D));
sqrtEigenvals(sqrtEigenvals < 0) = 0;

sqrtRL = P*sqrtEigenvals*P^(-1);

krausOpX = sqrt(pGenCondSingleandA(2))*zket(2,1)*zket(3,2)' + sqrt(pGenCondSingleandA(3))*zket(2,2)*zket(3,3)'; 
%krausOpX = kron(krausOpX, sqrtm(GammaR + GammaL)); %only comp basis used for key gen
krausOpX = kron(krausOpX, sqrtRL);
krausOps = {krausOpX}; 
newParams.krausOps = krausOps; 

krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausOps{index}'*krausOps{index};
end
debugInfo.storeInfo("krausSum", krausSum);

%% Z map 
proj0 = kron(diag([1,0]), eye(dimB) ); %ident on B' system (B - flag states)
proj1 = kron(diag([0,1]), eye(dimB) );
keyProj = {proj0,proj1};

%% set key map, kraus ops, and block diagonal structure in new parameters
newParams.krausOps = krausOps;
newParams.keyProj = keyProj;
newParams.blockDimsA = 3;
newParams.blockDimsB = [3, ones([1, 7])]; 

%%
newParams.n = floor(params.pGen*params.N); %gen
newParams.m = floor((1-params.pGen)*params.N); %test

newParams.coarseGrain = eye([dimA*dimB,dimA*dimB]);
newParams.Sigma = size(newParams.coarseGrain, 1);
newParams.Lambda = size(newParams.coarseGrain, 2); 

n_sift = floor(newParams.n*pxB); %%only used for optimal epsilons (does not influence security)
t = max(params.tau(:)); 

epsilons = optimalEpsVals(params.N, n_sift, t, dimA, params.epsSound, params.epsATfrac);
%epsilons = [params.epsSound/4,params.epsSound/4,params.epsSound/4,params.epsSound/4]; 
newParams.epsilons = epsilons; 
end


function mustBeProbDistCell(input)
mustBeProbDist([input{:}])
end

function mustBePOVM(S)
    arguments
        S (1, :) cell
    end
    ident = zeros(size(S{1})); 
    for i = 1:numel(S)
        ident = ident + S{i}; 
    end
    if ~all(ismembertol(real(ident), eye(size(ident))), 'all') || ~all(ismembertol(imag(ident), zeros(size(ident)),"Datascale",1),'all')
        throw(MException('mustBePOVM:NotPOVM', ...
            'The given cell aray does not sum to identity.'));
    end
end
