function qkdInput = ThreeStateDataTableDifInt_Preset()
% BasicBB84_LossyPreset a preset for a simple BB84 protocol with a
% lossy channel. Loss is included as a third dimension orthogonal to
% Alice's qubit dimnensions. Alice has been reduced to two dimensions from
% four via the Schmidt decomposition.
% Reviewed by Devashish Tupkary 2023/09/18.

qkdInput = QKDSolverInput();

%% Parameters

%Load data from file

% %Run 3
matdata = load("Data/ohioDataVarsObsTableTopRun3.mat");

% %Run 4
% matdata = load("Data/ohioDataVarsObsTableTopRun4.mat");

%Run 3 and 4
% matdata = load("Data/ohioDataVarsObsTableTopRun3and4.mat");

%add data to qkdInput
%Predicted probabilities
qkdInput.addFixedParameter("ohioData", matdata.ohioProbDist)
%Observation Order
qkdInput.addFixedParameter("ohioOrder", matdata.ohioOrder)
%G matrix
qkdInput.addFixedParameter("ohioG", matdata.ohioG)
%Observations
qkdInput.addFixedParameter("ohioObservations", matdata.ohioData) %used to compute key with
%Total signals sent
qkdInput.addFixedParameter("N", matdata.totalSignals);


%Add misalignment
% qkdInput.addFixedParameter("misalignmentAngle",0.1);
%Example of adding misalingment as a scan parameter
%qkdInput.addScanParameter("misalignmentAngle",num2cell(linspace(0,pi/4,11)));

%Add depolarization
% qkdInput.addFixedParameter("depolarization", 0.1);

%Add dark count rate
% qkdInput.addFixedParameter("darkCountRate", 1e-5);

%Add loss
lossdB = linspace(0, 10, 3);
lossEta = 10.^(-lossdB/10);
%qkdInput.addScanParameter("eta", num2cell(lossEta));
qkdInput.addFixedParameter('eta', 10^(-1.5));

%Define basis choice probability of Alice and Bob (Z basis)
qkdInput.addFixedParameter("pGen", 0.9);
%qkdInput.addScanParameter("pGen", num2cell(linspace(0.1, 0.9, 21))); 

%Only used in channel to generate click pattern order (any other value could be used as well) 
qkdInput.addFixedParameter("pzB", 1/2);

%Error correction efficiency f >= 1, f=1 is at Shanon limit 
qkdInput.addFixedParameter("f",1.16);

%% Finite Correction parameters 
qkdInput.addFixedParameter("epsSound", 1e-8);
qkdInput.addFixedParameter("epsATfrac", 0.95); 

t = matdata.t;
qkdInput.addFixedParameter("t", t); 
tau = matdata.tau; %initialize array for t values
qkdInput.addFixedParameter("tau", tau*t); 
tsift = matdata.tsift;
qkdInput.addFixedParameter("tsift", tsift);

%% Decoy Parameters

% Add loaded data to preset
qkdInput.addFixedParameter("GROUP_decoysSignal_1", matdata.intensities(3)); %signal intensity H
qkdInput.addFixedParameter("GROUP_decoysSignal_2", matdata.intensities(1)); %signal intensity R
qkdInput.addFixedParameter("GROUP_decoysSignal_3", matdata.intensities(2)); %signal intensity L

qkdInput.addFixedParameter("GROUP_decoys1_1", matdata.intensities(6));% decoy intensity 1 H
qkdInput.addFixedParameter("GROUP_decoys1_2", matdata.intensities(4));% R
qkdInput.addFixedParameter("GROUP_decoys1_3", matdata.intensities(5));% L

qkdInput.addFixedParameter("GROUP_decoys2_1", matdata.intensities(9));% decoy intensity 2 H (something slightly above 0 is usually optimal.)
qkdInput.addFixedParameter("GROUP_decoys2_2", matdata.intensities(7));% R
qkdInput.addFixedParameter("GROUP_decoys2_3", matdata.intensities(8));% L

%Decoy probabilities
qkdInput.addFixedParameter("GROUP_decoyProbs_1", matdata.decoyProbs(1)); %p(mu1|test)
qkdInput.addFixedParameter("GROUP_decoyProbs_2", matdata.decoyProbs(2));
qkdInput.addFixedParameter("GROUP_decoyProbs_3", matdata.decoyProbs(3));


%% Modules
% description 
descriptionModule = QKDDescriptionModule(@ThreeStateDataTableDifInt_Desc);
qkdInput.setDescriptionModule(descriptionModule);

% channel model
channelModule = QKDChannelModule(@ThreeStateDataDifInt_ChannelFunc); %Channel for decoys
qkdInput.setChannelModule(channelModule);

% key rate function
keyModuleOptions = struct(); 
keyModuleOptions.decoyTolerance = 1e-14;
keyModuleOptions.decoySolver = "Mosek";
keyModuleOptions.decoyForceSep = true;
keyModule = QKDKeyRateModule(@ThreeStateDataTableDifInt_KeyRateFunc, keyModuleOptions);
qkdInput.setKeyRateModule(keyModule);

% optimization
optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",1),struct("verboseLevel",1));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 20;
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.blockDiagonal = true;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

%% global options
%ErrorHandling.CatchWarn
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.DontCatch,"verboseLevel",1,"cvxSolver","sdpt3", "cvxPrecision", "default"));