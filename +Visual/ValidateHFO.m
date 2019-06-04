function [nEventValidity, strEventValidity] = ValidateHFO(rhfo, frhfo , ElectrodeInd)
%% %%%%%%%%%%%%%% This comes from ECE
%%Paths
strPaths.Main = '\\fl-daten\NCH_Forschungen\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Katja Rutz\Intraoperative HFO Validation HD 2\';
strPaths.HFOAnalysisResults = '\\fl-daten\NCH_Forschungen\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Maxine Schreiber\HFO Analysis Results\HFO Analysis 1 Extracted Data Added Channels 190111\' ; % [strPaths.Main,'HFO Analysis Results\'];
strPaths.Data = '\\fl-daten\NCH_Forschungen\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Maxine Schreiber\Results\Extracted Data Added Channels\'; % [strPaths.Main,'Data\'];
strPaths.HFOValidationResults = [strPaths.Main,'HFO Validation Results\'];
strPaths.ImagesExamplesOfEvents = [strPaths.Main,'Examples of Events\'];
addpath(genpath(strPaths.Main))

% File names for saving results
load('\\fl-daten\NCH_Forschungen\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Maxine Schreiber\code\ece\FIR_2KHz.mat')
addpath(genpath('\\fl-daten\NCH_Forschungen\NCH_FL_Forschungsprojekte\Epilepsy\_Tommaso\hfEEG\hfEEG software analysis\hfEEGAquistionRoutines\'));
%% %%%%%%%%%%%%%%%%%
nbSamples = rhfo.Data.nbSamples;
HFOARSI_ToPlot                           = getMarkingsFR1(frhfo, rhfo);
HFOARSI_ToPlot.MarkingsFR1_ToPlot        = HFOARSI_ToPlot.MarkingsFR1;
VInterfacePara.nChan                     = rhfo.Data.nbChannels;                                 %  scalar value number of channels
VInterfacePara.strEventType              = 'fast ripple';                           %  string 'ripple' or 'fast ripple'
VInterfacePara.ListOfChannelsToPlot      =  ElectrodeInd;                           %  indeces set by hand
VInterfacePara.HFOAnalysisResultsAllData = HFOARSI_ToPlot;
VInterfacePara.ElectrodeLabels           = strrep(rhfo.Data.channelNames','_','\_');%  cell of chars i.e. 'channel label'
VInterfacePara.y_shiftRaw                = 500;                                     %  shift in plot so that data is not overlapped
VInterfacePara.y_shiftFiltered           = 10;                                      %  shift in plot so that data is not overlapped
VInterfacePara.strSaveImagesFolderName   = 'E:\GENEVA\Processed Data\Geneva\';      %  string of save path
VInterfacePara.t_axis                    = 1:nbSamples;                             %  1:numberOfSamples
%%
% [nEventValidity, strEventValidity]       = Event_Validation_GUI_190304(VInterfacePara);
[nEventValidity, strEventValidity]       = Event_Validation_GUI_190304(ValidationInterfaceParams);
end

function HFOARSI_ToPlot = getMarkingsFR1(rhfo, frhfo)
    nbChan    = rhfo.Data.nbChannels;
    fs        = rhfo.Data.sampFreq; 
    
    HFOARSI_ToPlot.MarkingsFR1 = [];
    marks = rhfo.Events.Markings;
    for iChan = 1:nbChan
        HFOARSI_ToPlot.HFOAnalysisResultsAllChannels{iChan,1}.signal          = rhfo.Data.signal(:,iChan)';
        HFOARSI_ToPlot.HFOAnalysisResultsAllChannels{iChan,1}.signalFilt      = rhfo.filtSig.filtSignal(:,iChan)';
        HFOARSI_ToPlot.HFOAnalysisResultsAllChannels{iChan,1}.signalFiltFR    = frhfo.filtSig.filtSignal(:,iChan)';
        nbEvents = length(marks.start{iChan});

        EventBlock = ([fs*ones(nbEvents,1)*iChan, marks.start{iChan}', marks.end{iChan}', marks.len{iChan}']./fs);
        HFOARSI_ToPlot.MarkingsFR1  = [HFOARSI_ToPlot.MarkingsFR1 ; EventBlock];
    end
end
