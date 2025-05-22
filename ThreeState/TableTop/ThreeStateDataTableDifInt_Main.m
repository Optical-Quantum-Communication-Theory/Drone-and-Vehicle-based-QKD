clear all

%pick preset
qkdInput = ThreeStateDataTableDifInt_Preset();
%inside preset select data of experimental run

%run the QKDSolver with this input and store results
results = MainIteration(qkdInput);

%save the results and preset to a file.
% save("data.mat","results","qkdInput");

%% plot the result
%QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")