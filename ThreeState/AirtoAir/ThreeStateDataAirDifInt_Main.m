clear all

%% pick preset

%Flight 1 Nov 13
% qkdInput = ThreeStateDataAirDifInt_Preset_Flight1Nov13();

%Flight 1 Nov 2nd
qkdInput = ThreeStateDataAirDifInt_Preset_Flight1Nov2();

%Flight 2 Nov 2nd
% qkdInput = ThreeStateDataAirDifInt_Preset_Flight2Nov2();


%Windowed Flight 1 Nov 2nd
% qkdInput = ThreeStateDataAirDifInt_Preset_WindowedFlight1Nov2();

%Windowed Flight 2 Nov 2nd
% qkdInput = ThreeStateDataAirDifInt_Preset_WindowedFlight2Nov2();

%Flight 1 and 2 combined
% qkdInput = ThreeStateDataAirDifInt_Preset_CombinedNov2();


%% run the QKDSolver with this input and store results
results = MainIteration(qkdInput);

%save the results and preset to a file.
% save("data.mat","results","qkdInput");

%% plot the result
%QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")