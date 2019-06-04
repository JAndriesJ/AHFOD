% This script is create for the processing of a sub selection of the data
% uploaded in the SciRep Paper.
clear
%% Run Detector for Ripples
warning('off','all')
%% Load parameters and data
% Patient numbers in the publication
% patN_i  =  [7, 17,  9,  18,  4, 10,  1  14  11   2  12  13   8  16  15   3  19  20   5   6];
% Patient numbers as on server
% patN_ece = [1,  3,  5,   7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23];

PatNum = {'01' '05' '11' '12' '13' '16' '22' '23'};% as on the server
for iFile = 1:8
    disp(['Currently running patient: ',num2str(PatNum{iFile})])
    RParaPath = [pwd, '\+Demo\ZurichDemo1Morph\Parameters\RMorphPara.mat'];
    FRParaPath = [pwd, '\+Demo\ZurichDemo1Morph\Parameters\FRMorphPara.mat'];
    DataDir  = ['E:\GENEVA\Processed Data\Zurich\Data selection for sciRep\Pat',num2str(PatNum{iFile}),'\'];
    %% Run the detecor
    hfodet = Detections;
    hfodet = hfodet.setPDPaths(RParaPath, FRParaPath, DataDir);
    hfodet = hfodet.runDetector('morph', false, 3 , '', true, 'RipAndFRip');
    hfodet.exportHFOsummary
end
%% Sanity check
% ReffPath = 'E:\GENEVA\Processed Data\Zurich\Old Results selection for sciRep\Patient 05 FC\HFO_Analysis_Results_Patient_05_Recording_01_Interval_101.mat';
% load(ReffPath,'HFOAnalysisResults');
% OutputPath = 'E:\GENEVA\Processed Data\Zurich\Data selection for sciRep\Pat05\RippleHFOSummary\HFOSummaryMat.mat'; %Ripple
% OutputPath = 'E:\GENEVA\Processed Data\Zurich\Data selection for sciRep\Pat05\FastRippleHFOSummary\HFOSummaryMat.mat'; %Fast Ripple
% load(OutputPath);
%
% clc
% chan = 9;
% % ReferenceResult = HFOAnalysisResults.HFOAnalysisResultsAllChannels{chan, 1}.MarkingsRipple1(:,2)';
% ReferenceResult = HFOAnalysisResults.HFOAnalysisResultsAllChannels{chan, 1}.MarkingsFR1(:,3)';
% ComputedResult = HFOSummaryMat.hfoInfo{1, 1}.Events.Markings.end{1, chan}/2000 ; 
%%

