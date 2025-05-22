function qkdInput = ThreeStateDataAirtoCarDifInt_Preset_Run2()
% BasicBB84_LossyPreset a preset for a simple BB84 protocol with a
% lossy channel. Loss is included as a third dimension orthogonal to
% Alice's qubit dimnensions. Alice has been reduced to two dimensions from
% four via the Schmidt decomposition.
% Reviewed by Devashish Tupkary 2023/09/18.

qkdInput = QKDSolverInput();

%% Parameters

% Run 2 Nov 14
matdata = load("Data/ohioDataVarsObsAirtoCarRun2Nov14.mat");

%add data to qkdInput
qkdInput.addFixedParameter("ohioData", matdata.ohioProbDistCar) %Probability distribution only used to compare data with
qkdInput.addFixedParameter("ohioOrder", matdata.ohioOrderCar)
qkdInput.addFixedParameter("ohioG", matdata.ohioGCar)
qkdInput.addFixedParameter("ohioObservations", matdata.ohioDataCar) %used to compute key with

%Add misalignment
% qkdInput.addFixedParameter("misalignmentAngle",0.1);
%Example of adding misalingment as a scan parameter
%qkdInput.addScanParameter("misalignmentAngle",num2cell(linspace(0,pi/4,11)));

%Add depolarization
% qkdInput.addFixedParameter("depolarization", 0.1);

%Add dark count rate
% qkdInput.addFixedParameter("darkCountRate", 1e-5);

%Add loss
% lossdB = linspace(0, 10, 3);
% lossEta = 10.^(-lossdB/10);
%qkdInput.addScanParameter("eta", num2cell(lossEta));
% qkdInput.addFixedParameter('eta', 10^(-1.5));

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

% %Run 2 Nov 14
qkdInput.addFixedParameter("N", 1284704352);

% t parameter
t = 10^(-10);
qkdInput.addFixedParameter("t", t); 
tau = ones([3, 8, 3]); %initialize array for t values

% % %Run 2 Nov 14
tau(2,2,3) = 1e4;
tau(3,2,3) = 1e4;
tau(3,3,3) = 1e4;

qkdInput.addFixedParameter("tau", tau*t); 
qkdInput.addFixedParameter("tsift", 0) 

%% Decoy Parameters

%Run 2 Nov 14
qkdInput.addFixedParameter("GROUP_decoysSignal_1", 0.44); %signal intensity H
qkdInput.addFixedParameter("GROUP_decoysSignal_2", 0.44); %signal intensity R
qkdInput.addFixedParameter("GROUP_decoysSignal_3", 0.44); %signal intensity L

qkdInput.addFixedParameter("GROUP_decoys1_1", 0.16383); % decoy intensity 1 H
qkdInput.addFixedParameter("GROUP_decoys1_2", 0.18263); % R
qkdInput.addFixedParameter("GROUP_decoys1_3", 0.17505);% L

qkdInput.addFixedParameter("GROUP_decoys2_1", 0.0079147); % decoy intensity 2 H 
qkdInput.addFixedParameter("GROUP_decoys2_2", 0.010329); % R
qkdInput.addFixedParameter("GROUP_decoys2_3", 0.014392); % L

%Decoy probs are equal for both runs
qkdInput.addFixedParameter("GROUP_decoyProbs_1", 4821/9419); %p(mu1|test)
qkdInput.addFixedParameter("GROUP_decoyProbs_2", 5929/18838);
qkdInput.addFixedParameter("GROUP_decoyProbs_3", 3267/18838);

%% Modules
% description 
descriptionModule = QKDDescriptionModule(@ThreeStateDataAirtoCarDifInt_Desc);
qkdInput.setDescriptionModule(descriptionModule);

% channel model
channelModule = QKDChannelModule(@ThreeStateDataDifInt_ChannelFunc); %Channel for decoys
qkdInput.setChannelModule(channelModule);

% key rate function
keyModuleOptions = struct(); 
keyModuleOptions.decoyTolerance = 1e-14;
keyModuleOptions.decoySolver = "Mosek";
keyModuleOptions.decoyForceSep = true;
keyModule = QKDKeyRateModule(@ThreeStateDataAirtoCarDifInt_KeyRateFunc, keyModuleOptions);
qkdInput.setKeyRateModule(keyModule);

% optimization
optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",1),struct("verboseLevel",1));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 25;
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.blockDiagonal = true;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

%% global options
%ErrorHandling.CatchWarn
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.DontCatch,"verboseLevel",1,"cvxSolver","mosek", "cvxPrecision", "default"));