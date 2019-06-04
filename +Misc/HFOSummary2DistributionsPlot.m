close
clear
clc
%%
DirID = {'1757\224081' '1785\236222' '1785\242175' '1786\221485' '1786\225402' '1806\231685' '1806\231697' '1823\254172' '1823\254169' '2012\255619' '2012\257707' '2091\247147' '2091\247134'};
for iDir = 13:13
saveDir = ['E:\GENEVA\Processed Data\Geneva\PAT_',DirID{iDir},'\Data\HFOSummary\'];
% saveDir = ['E:\GENEVA\Processed Data\Geneva\PAT_',DirID{iDir}(1:4),'\'];
dataDir = [saveDir, 'HFOSummaryMat.mat'];
load(dataDir)

nbFiles   = length(HFOSummaryMat.Ripples);
nbChan    = HFOSummaryMat.Ripples{1}.Data.nbChannels;
chanNames = HFOSummaryMat.Ripples{1}.Data.channelNames;

[RHFOraMat, FRHFOraMat, RandFRHFOraMat] = getRateMats(HFOSummaryMat, nbFiles, nbChan);

%%% get masks
maskNonLateralContacts = ~contains(chanNames, {'10','11','12','13','14','15','16','17','18','19','20'} );
S_RHFOraMat      = RHFOraMat(:,maskNonLateralContacts); 
S_FRHFOraMat     = FRHFOraMat(:,maskNonLateralContacts); 
S_RandFRHFOraMat = RandFRHFOraMat(:,maskNonLateralContacts);
S_chanNames      = chanNames(maskNonLateralContacts);


[HFOArea, maskHFOArea,  valMaxTHR] = Detections.GetHFOAreaMat(S_RandFRHFOraMat');
%%% saveSummary(S_RHFOraMat, S_FRHFOraMat, S_RandFRHFOraMat, HFOArea, valMaxTHR, dataDir);
Detections.exportHFOrateDistribution(S_RHFOraMat, S_FRHFOraMat, S_RandFRHFOraMat, HFOArea, valMaxTHR, S_chanNames, saveDir, 'noLat')
Detections.exportHFOrateDistribution( [], [], S_RandFRHFOraMat, HFOArea, valMaxTHR, S_chanNames, saveDir, 'R&FRnoLat')
Detections.exportReproducibilityPlots(S_RandFRHFOraMat, HFOArea, maskHFOArea, S_chanNames, saveDir);
%% %%%
[RNoiseMat, RBaseLMat, FRNoiseMat, FRBaseLMat] = getNoiseMats(HFOSummaryMat, nbFiles, nbChan);
S_RNoiseMat  = RNoiseMat(:,maskNonLateralContacts);
S_RBaseLMat  = RBaseLMat(:,maskNonLateralContacts);
S_FRNoiseMat = FRNoiseMat(:,maskNonLateralContacts);
S_FRBaseLMat = FRBaseLMat(:,maskNonLateralContacts);

Detections.exportNoiseDist( S_RBaseLMat,  S_FRBaseLMat,  S_RNoiseMat,  S_FRNoiseMat, HFOArea,  S_chanNames, saveDir)
end
%%

function [RHFOraMat, FRHFOraMat, RandFRHFOraMat] = getRateMats(HFOSummaryMat, nbFiles, nbChan)
    RHFOraMat      = zeros(nbFiles,nbChan);
    FRHFOraMat     = zeros(nbFiles,nbChan);
    RandFRHFOraMat = zeros(nbFiles,nbChan);

    for iFile = 1:nbFiles
        RHFOraMat(iFile,:)      = HFOSummaryMat.Ripples{iFile}.Events.Rates;
        FRHFOraMat(iFile,:)     = HFOSummaryMat.FastRipples{iFile}.Events.Rates;
        RandFRHFOraMat(iFile,:) = HFOSummaryMat.RandFR{iFile}.Rates.RippleANDFastRipple';
    end
end

function [RNoiseMat, RBaseLMat, FRNoiseMat, FRBaseLMat] = getNoiseMats(HFOSummaryMat, nbFiles, nbChan)
    RNoiseMat  = zeros(nbFiles,nbChan);
    RBaseLMat  = zeros(nbFiles,nbChan);
    FRNoiseMat = zeros(nbFiles,nbChan);
    FRBaseLMat = zeros(nbFiles,nbChan);
    for iFile = 1:nbFiles
        RNoiseMat(iFile,:)    = HFOSummaryMat.Ripples{1, iFile}.baseline.maxNoisemuV;
        RBaseLMat(iFile,:)    = HFOSummaryMat.Ripples{1, iFile}.baseline.baselineThr ;
        FRNoiseMat(iFile,:)   = HFOSummaryMat.FastRipples{1, iFile}.baseline.maxNoisemuV;
        FRBaseLMat(iFile,:)   = HFOSummaryMat.FastRipples{1, iFile}.baseline.baselineThr;
    end
end

function [] = saveSummary(RHFOraMat, FRHFOraMat, RandFRHFOraMat, HFOArea, valMaxTHR, dataDir)
HFOSummaryMat.Rates.RHFOMat      = RHFOraMat;
HFOSummaryMat.Rates.FRHFOMat     = FRHFOraMat;
HFOSummaryMat.Rates.RandFRHFOMat = RandFRHFOraMat;
HFOSummaryMat.Area.HFOArea          = HFOArea;
HFOSummaryMat.Area.HFOthreshold     = valMaxTHR;
save(dataDir ,'HFOSummaryMat')
end
