clear all
%pick preset
qkdInput = ThreeStateDataGroundDifInt_Preset();

%Load data
% Run 1 Nov 14
matdata = load("Data/ohioDataVarsObsGroundRun1.mat");

%add data to qkdInput
qkdInput.addFixedParameter("ohioData", matdata.ohioProbDist) %Probability distribution only used to compare data with
qkdInput.addFixedParameter("ohioOrder", matdata.ohioOrder)
qkdInput.addFixedParameter("ohioG", matdata.ohioG)
qkdInput.addFixedParameter("ohioObservations", matdata.ohioData) %used to compute key with

%run the QKDSolver with this input and store results
results = MainIteration(qkdInput);

%save the results and preset to a file.
% save("data.mat","results","qkdInput");

%% plot the result
%QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")