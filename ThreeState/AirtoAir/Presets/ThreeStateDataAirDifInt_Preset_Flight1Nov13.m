function qkdInput = ThreeStateDataAirDifInt_Preset_Flight1Nov13()
% BasicBB84_LossyPreset a preset for a simple BB84 protocol with a
% lossy channel. Loss is included as a third dimension orthogonal to
% Alice's qubit dimnensions. Alice has been reduced to two dimensions from
% four via the Schmidt decomposition.
% Reviewed by Devashish Tupkary 2023/09/18.

qkdInput = QKDSolverInput();

%% Parameters

%Flight 1 Nov 13
matdata = load("Data/ohioDataVarsObsAirFlight1Nov13.mat");

%add data to qkdInput
qkdInput.addFixedParameter("ohioData", matdata.ohioProbDistAir)
qkdInput.addFixedParameter("ohioOrder", matdata.ohioOrderAir)
qkdInput.addFixedParameter("ohioG", matdata.ohioGAir)
qkdInput.addFixedParameter("ohioObservations", matdata.ohioDataAir) %used to compute key with

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

% Flight 1 Nov 13
qkdInput.addFixedParameter("N", 1081348008);

%Flight 1 Nov 13
t = 10^(-8);
qkdInput.addFixedParameter("t", t); 
tau = ones([3, 8, 3]); %initialize array for t values
tau(1,2,3) = 1e5;
tau(2,2,3) = 1e4;
tau(2,4,3) = 1e4;
qkdInput.addFixedParameter("tau", tau*t); 
qkdInput.addFixedParameter("tsift", 0);

%% Decoy Parameters

%Flight 1 Nov 13
qkdInput.addFixedParameter("GROUP_decoysSignal_1", 0.18); %signal intensity H
qkdInput.addFixedParameter("GROUP_decoysSignal_2", 0.18); %signal intensity R
qkdInput.addFixedParameter("GROUP_decoysSignal_3", 0.18); %signal intensityL

qkdInput.addFixedParameter("GROUP_decoys1_1", 0.079526);% decoy intensity 1 H
qkdInput.addFixedParameter("GROUP_decoys1_2", 0.070573);% R
qkdInput.addFixedParameter("GROUP_decoys1_3", 0.077108);% L

qkdInput.addFixedParameter("GROUP_decoys2_1", 0.012227);% decoy intensity 2 H (something slightly above 0 is usually optimal.)
qkdInput.addFixedParameter("GROUP_decoys2_2", 0.0060329);% R
qkdInput.addFixedParameter("GROUP_decoys2_3", 0.01191);% L

%Decoy probabilities
qkdInput.addFixedParameter("GROUP_decoyProbs_1", 4821/9419); %p(mu1|test)
qkdInput.addFixedParameter("GROUP_decoyProbs_2", 5929/18838);
qkdInput.addFixedParameter("GROUP_decoyProbs_3", 3267/18838);

%% Modules
% description 
descriptionModule = QKDDescriptionModule(@ThreeStateDataAirDifInt_Desc_FlightNov13);
qkdInput.setDescriptionModule(descriptionModule);

% channel model
channelModule = QKDChannelModule(@ThreeStateDataDifInt_ChannelFunc); %Channel for decoys
qkdInput.setChannelModule(channelModule);

% key rate function
keyModuleOptions = struct(); 
keyModuleOptions.decoyTolerance = 1e-14;
keyModuleOptions.decoySolver = "Mosek";
keyModuleOptions.decoyForceSep = true;
keyModule = QKDKeyRateModule(@ThreeStateDataAirDifInt_KeyRateFunc, keyModuleOptions);
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
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.DontCatch,"verboseLevel",1,"cvxSolver","mosek", "cvxPrecision", "default"));