clear all

%% pick preset

% Run 1 
% qkdInput = ThreeStateDataAirtoCarDifInt_Preset_Run1();

%Run 2 Nov 14
qkdInput = ThreeStateDataAirtoCarDifInt_Preset_Run2();

%% run the QKDSolver with this input and store results
results = MainIteration(qkdInput);

%save the results and preset to a file.
% save("data.mat","results","qkdInput");

%% plot the result
%QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")