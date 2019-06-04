%% Close all figures, clear variables and command window
close all
% clear all
% clc

% Line 34 to change file number for each run
% Line 46 to change event type (1:Ripple, 2:FR)
% Ctrl+C to stop running without saving the results
% From line 129 change order of channels to show
%% Paths
strPaths.Main = '\\fl-daten\NCH_Forschungen\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Katja Rutz\Intraoperative HFO Validation HD 2\';
strPaths.HFOAnalysisResults = '\\fl-daten\NCH_Forschungen\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Maxine Schreiber\HFO Analysis Results\HFO Analysis 1 Extracted Data Added Channels 190111\' ; % [strPaths.Main,'HFO Analysis Results\'];
strPaths.Data = '\\fl-daten\NCH_Forschungen\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Maxine Schreiber\Results\Extracted Data Added Channels\'; % [strPaths.Main,'Data\'];
strPaths.HFOValidationResults = [strPaths.Main,'HFO Validation Results\'];
strPaths.ImagesExamplesOfEvents = [strPaths.Main,'Examples of Events\'];

% Add all subfolders to path
addpath(genpath(strPaths.Main))

% File names for saving results
strFormat.SaveValidationResults = '%s_Validation_Results_File_%s_Group_%.2d_Channel_%.2d_%s.mat';
strFormat.SaveSnapshots = strPaths.ImagesExamplesOfEvents;
strFormat.ValidationComplete = '%s Validation complete for Patient %s Recording %d on Day %d for Channel %d %s';

load('\\fl-daten\NCH_Forschungen\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Maxine Schreiber\code\ece\FIR_2KHz.mat')

results_event_table_validated_dir = '\\fl-daten\NCH_Forschungen\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Maxine Schreiber\Results Event Table Validated Stage 2 180927\';
addpath(genpath('\\fl-daten\NCH_Forschungen\NCH_FL_Forschungsprojekte\Epilepsy\_Tommaso\hfEEG\hfEEG software analysis\hfEEGAquistionRoutines\'));

%% List of files in folder 'HFO Analysis Results'
strListOfFiles = dir([strPaths.HFOAnalysisResults,'*.mat']);
strListOfFiles = {strListOfFiles.name}';
nNumberOfFiles = size(strListOfFiles,1);

%% Load each file
% for nFile = 3:nNumberOfFiles
for nFile = 4
    % size(strListOfFiles,1) %% CHANGE THE NUMNER OF FILE YOU ARE ANALYZING BEFORE RUNNING THE SCRIPT
    % Choose event type
    for nEventType = 2 % CHANGE THIS NUMBER TO SELECT EVENT TYPE (1: Ripple, 2: Fast ripple)
        % Select group of electrodes
        nGroup = 5;
        %%
        strFileName = strListOfFiles{nFile};
        FolderName = strFileName(1:end-4);
        strValidationResultsFolder = [strPaths.HFOValidationResults,strFileName(1:end-4),'\'];
        mkdir(strValidationResultsFolder)
        strFormat.SaveSnapshots = [strPaths.ImagesExamplesOfEvents,strFileName(1:end-4),'\'];
        mkdir(strFormat.SaveSnapshots)
        strFileFullPath = [strPaths.HFOAnalysisResults,strFileName];
        Vars_HFO = load(strFileFullPath);
        dataAll = load([strPaths.Data,strrep(strFileName,'HFO_','')]);
        dataAll = dataAll.data;
        t_axis = dataAll.t;
        
        %% Select data from the group
        indGroupElectrodes = dataAll.nElectrodeGroups==nGroup;
        nGroupElectrodes = find(indGroupElectrodes);
        nNewElectrodes = nGroupElectrodes-nGroupElectrodes(1)+1;
        %% Coordinates
        coord_temp = dataAll.Coord_new(nGroupElectrodes,:);
        clear Datasetup Dist
        for i = 1:size(coord_temp,1)
            for j = 1:size(coord_temp,1)
                Dist(i,j) = norm(coord_temp(i,:)-coord_temp(j,:));
            end
            [Datasetup(i).Dist_val, Datasetup(i).Dist_ord] = sort(Dist(i,:));
        end
        %%
        switch nEventType
            case 1
                HFO_FR = Vars_HFO.HFO_R;
            case 2
                HFO_FR = Vars_HFO.HFO_FR;
        end
        HFOobj = HFO_FR.HFOobj;
        %% Get data for selected group
        EVENTS3 = HFOobj.EVENTS3;
        EVENTS3 = EVENTS3(ismember(EVENTS3(:,1),nGroupElectrodes),:);
        % Replace electrode numbers
        for nElectrode = 1:length(nGroupElectrodes)
            EVENTS3(EVENTS3(:,1)==nGroupElectrodes(nElectrode),1) = -nElectrode;
        end
        EVENTS3(:,1) = -EVENTS3(:,1);
        
        %%
        Markings = EVENTS3(:,1:3);
        Markings(:,2) = (Markings(:,2)-1)/2000;
        Markings(:,3) = (Markings(:,3)-1)/2000;
        Markings(:,4) = Markings(:,2)+Markings(:,3);
        Markings = Markings(:,[1,2,4,3]);
        %% Order in time
        [~,indSort] = sort(Markings(:,2),'ascend');
        Markings = Markings(indSort,:);
        %%
        ElectrodeLabels = dataAll.ElectrodeLabels_new(nGroupElectrodes);
        nNumberOfChannels = length(ElectrodeLabels);
        
        %%
        HFOAnalysisResultsAllChannels = cell(nNumberOfChannels,1);
        for nChannel = 1:nNumberOfChannels
            HFOAnalysisResultsAllChannels{nChannel}.signal = dataAll.x_bip_new(nGroupElectrodes(nChannel),:);
            HFOAnalysisResultsAllChannels{nChannel}.signal = HFOAnalysisResultsAllChannels{nChannel}.signal...
                -mean(HFOAnalysisResultsAllChannels{nChannel}.signal);
            b = filter.Rb;
            a = filter.Ra;
            HFOAnalysisResultsAllChannels{nChannel}.signalFilt = filtfilt(b,a,dataAll.x_bip_new(nGroupElectrodes(nChannel),:));
            b = filter.FRb;
            a = filter.FRa;
            HFOAnalysisResultsAllChannels{nChannel}.signalFiltFR = filtfilt(b,a,dataAll.x_bip_new(nGroupElectrodes(nChannel),:));
        end
        
        
        nChannelRange_ToDelete = [];
        switch strFileName
            case 'HFO_pat161031_CD_edf8.mat'
                switch nGroup
                    case 3
                        nChannelRange_ToDelete = 4;
                    case 6
                        nChannelRange_ToDelete = [3,4,6];
                end
            case 'HFO_pat170515_CD_edf5.mat'
                switch nGroup
                    case 3
                        nChannelRange_ToDelete = 4;
                end
            case 'HFO_pat170710_CD_edf4.mat'
                switch nGroup
                    case 6
                        nChannelRange_ToDelete = [4,5,6];
                end
            case 'HFO_pat171024_CD_edf1_POST1_1.mat'
                switch nGroup
                    case 4
                        nChannelRange_ToDelete = 5;
                end
            case 'HFO_pat171024_CD_edf4_POST3_11.mat'
                switch nGroup
                    case 3
                        nChannelRange_ToDelete = [5,6];
                    case 4
                        nChannelRange_ToDelete = [3,5,6];
                    case 5
                        nChannelRange_ToDelete = [4,5,6];
                    case 6
                        nChannelRange_ToDelete = [3,5,6];
                end
            case 'HFO_pat180410_CD_edf1.mat'
                switch nGroup
                    case 3
                        nChannelRange_ToDelete = [1,2,3];
                    case 4
                        nChannelRange_ToDelete = [1,2];
                    case 6
                        nChannelRange_ToDelete = [1,2,4];
                end
            case 'HFO_pat160711_CD_edf3.mat' % Control+T
                switch nGroup
                    case 3
                        nChannelRange_ToDelete = [1,2,4];
                    case 4
                        nChannelRange_ToDelete = [1,2,3,4];
                    case 5
                        nChannelRange_ToDelete = [1,2];
                end
                
              
                %                 case 'xxx.mat' % Control+T
                %                 switch nGroup
                %                     case xxx
                %                         nChannelRange_ToDelete = [xxx];
                %                 end
                
        end
        for iChannel = 1:length(nChannelRange_ToDelete)
            HFOAnalysisResultsAllChannels{nChannelRange_ToDelete(iChannel)}.signal = 0*HFOAnalysisResultsAllChannels{nChannelRange_ToDelete(iChannel)}.signal;
            HFOAnalysisResultsAllChannels{nChannelRange_ToDelete(iChannel)}.signalFilt = 0*HFOAnalysisResultsAllChannels{nChannelRange_ToDelete(iChannel)}.signalFilt;
            HFOAnalysisResultsAllChannels{nChannelRange_ToDelete(iChannel)}.signalFiltFR = 0*HFOAnalysisResultsAllChannels{nChannelRange_ToDelete(iChannel)}.signalFiltFR;
        end
        HFOAnalysisResultsSingleInterval.HFOAnalysisResultsAllChannels = HFOAnalysisResultsAllChannels;
        
        %%
        ElectrodeLabels
        %% Markings to remove
        Markings_ToAdd = [];
        nChannelRange_ToDelete2 = [];
        switch strFileName
            case 'HFO_pat170515_CD_edf5.mat'
                switch nGroup
                    case 3
                        nChannelRange_ToDelete2 = 4;
                        Markings_ToAdd = Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:);
                        Markings_ToAdd = [Markings_ToAdd,7*ones(size(Markings_ToAdd,1),1)];
                        Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:) = [];
                end
            case 'HFO_pat170710_CD_edf4.mat'
                switch nGroup
                    case 6
                        nChannelRange_ToDelete2 = [4,5,6];
                        Markings_ToAdd = Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:);
                        Markings_ToAdd = [Markings_ToAdd,7*ones(size(Markings_ToAdd,1),1)];
                        Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:) = [];
                end
            case 'HFO_pat171024_CD_edf4_POST3_11.mat'
                switch nGroup
                    case 3
                        nChannelRange_ToDelete2 = [5,6];
                        Markings_ToAdd = Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:);
                        Markings_ToAdd = [Markings_ToAdd,7*ones(size(Markings_ToAdd,1),1)];
                        Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:) = [];
                    case 4
                        nChannelRange_ToDelete2 = [3,5,6];
                        Markings_ToAdd = Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:);
                        Markings_ToAdd = [Markings_ToAdd,7*ones(size(Markings_ToAdd,1),1)];
                        Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:) = [];
                    case 5
                        nChannelRange_ToDelete2 = [4,5,6];
                        Markings_ToAdd = Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:);
                        Markings_ToAdd = [Markings_ToAdd,7*ones(size(Markings_ToAdd,1),1)];
                        Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:) = [];
                    case 6
                        nChannelRange_ToDelete2 = [3,4,5,6];
                        Markings_ToAdd = Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:);
                        Markings_ToAdd = [Markings_ToAdd,7*ones(size(Markings_ToAdd,1),1)];
                        Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:) = [];
                end
                case 'HFO_pat180410_CD_edf1.mat'
                switch nGroup
                    case 6
                        nChannelRange_ToDelete2 = [1,2,4];
                        Markings_ToAdd = Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:);
                        Markings_ToAdd = [Markings_ToAdd,7*ones(size(Markings_ToAdd,1),1)];
                        Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:) = [];
                end
                
             case 'HFO_pat160711_CD_edf3.mat'
                switch nGroup
                    case 3
                        nChannelRange_ToDelete2 = [1,2,4];
                        Markings_ToAdd = Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:);
                        Markings_ToAdd = [Markings_ToAdd,7*ones(size(Markings_ToAdd,1),1)];
                        Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:) = [];
                    case 4
                        nChannelRange_ToDelete2 = [1,2,3,4];
                        Markings_ToAdd = Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:);
                        Markings_ToAdd = [Markings_ToAdd,7*ones(size(Markings_ToAdd,1),1)];
                        Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:) = [];
                    case 5
                        nChannelRange_ToDelete2 = [1,2];
                        Markings_ToAdd = Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:);
                        Markings_ToAdd = [Markings_ToAdd,7*ones(size(Markings_ToAdd,1),1)];
                        Markings(ismember(Markings(:,1),nChannelRange_ToDelete2),:) = [];
                end 
        end
        HFOAnalysisResultsSingleInterval.MarkingsFR1 = Markings;
        HFOAnalysisResultsSingleInterval.MarkingsFR1_ToPlot = Markings;
        
        %%
        switch nEventType
            case 1
                strEventType = 'ripple';
                y_shiftRaw = 750;
                y_shiftFiltered = 20;
                strEventType2 = 'Ripple';
                condSkip = isempty(Markings);
            case 2
                strEventType = 'fast ripple';
                strEventType2 = 'Fast_Ripple';
                y_shiftRaw = 750;
                y_shiftFiltered = 10;
                condSkip = isempty(Markings);
        end
        
        %% List of electrodes to display
        ListOfElectrodesForDisplay = [];
        if(length(ElectrodeLabels)>=5)
            for ch2plot = 1:length(ElectrodeLabels)
                ListOfElectrodesForDisplay = [ListOfElectrodesForDisplay;
                    ch2plot Datasetup(ch2plot).Dist_ord(2:5) Datasetup(ch2plot).Dist_ord(end)];
            end
        elseif(length(ElectrodeLabels)==4)
            for ch2plot = 1:length(ElectrodeLabels)
                ListOfElectrodesForDisplay = [ListOfElectrodesForDisplay;
                    ch2plot Datasetup(ch2plot).Dist_ord(2:4)];
            end
        elseif(length(ElectrodeLabels)==3)
            for ch2plot = 1:length(ElectrodeLabels)
                ListOfElectrodesForDisplay = [ListOfElectrodesForDisplay;
                    ch2plot Datasetup(ch2plot).Dist_ord(2:3)];
            end
        elseif(length(ElectrodeLabels)==2)
            for ch2plot = 1:length(ElectrodeLabels)
                ListOfElectrodesForDisplay = [ListOfElectrodesForDisplay;
                    ch2plot Datasetup(ch2plot).Dist_ord(2)];
            end
        elseif(length(ElectrodeLabels)==1)
            ListOfElectrodesForDisplay = 1;
        end
        
        %% Change some of the electrodes to display
        %         ListOfElectrodesForDisplay(1,:) = [1,2,4,3,5,12];
        %         ListOfElectrodesForDisplay(6,:) = [6,5,3,9,2,10];
        %         ListOfElectrodesForDisplay(7,:) = [7,8,10,4,11,3];
        %         ListOfElectrodesForDisplay(12,:) = [12,11,9,10,8,1];
        %%
        
        if(~condSkip)
            
            %% Loop over channels
            for nChan = 1 % 1:nNumberOfChannels
                
                ListOfChannelsToPlot = ListOfElectrodesForDisplay(nChan,:);
                %% Skip channel?
                if(~isempty(find(ismember(HFOAnalysisResultsSingleInterval.MarkingsFR1(:,1),ListOfChannelsToPlot),1)))
                    HFOAnalysisResultsSingleInterval_ToPlot = HFOAnalysisResultsSingleInterval;
                    HFOAnalysisResultsSingleInterval_ToPlot.MarkingsFR1 = ...
                        HFOAnalysisResultsSingleInterval_ToPlot.MarkingsFR1(ismember(HFOAnalysisResultsSingleInterval_ToPlot.MarkingsFR1(:,1),ListOfChannelsToPlot),:);
                    HFOAnalysisResultsSingleInterval_ToPlot.MarkingsFR1_ToPlot = HFOAnalysisResultsSingleInterval_ToPlot.MarkingsFR1;
                    
                    %% Run validation interface
                    ValidationInterfaceParams.nChan = nChan;
                    ValidationInterfaceParams.strEventType = strEventType;
                    ValidationInterfaceParams.ListOfChannelsToPlot = ListOfChannelsToPlot;
                    ValidationInterfaceParams.HFOAnalysisResultsAllData = HFOAnalysisResultsSingleInterval_ToPlot;
                    ValidationInterfaceParams.ElectrodeLabels = strrep(ElectrodeLabels,'_','\_');
                    ValidationInterfaceParams.y_shiftRaw = 500; %  y_shiftRaw;
                    ValidationInterfaceParams.y_shiftFiltered = 10;
                    ValidationInterfaceParams.strSaveImagesFolderName = strFormat.SaveSnapshots;
                    ValidationInterfaceParams.t_axis = t_axis;
                    %                                         error('Stop here')
                    [nEventValidity, strEventValidity] = ...
                        Event_Validation_GUI_190304(ValidationInterfaceParams);
                    
                    %% Exchange invalid with valid, artifact with invalid
                    %                     ind = find(nEventValidity==7);
                    %                     nEventValidity(ind) = 1;
                    %                     ind = find(nEventValidity==1);
                    %                     for ii = 1:length(ind)
                    %                         strEventValidity{ind(ii)} = 'Fast Ripple';
                    %                     end
                    %
                    %                     ind = find(nEventValidity==6);
                    %                     nEventValidity(ind) = 7;
                    %                     ind = find(nEventValidity==7);
                    %                     for ii = 1:length(ind)
                    %                         strEventValidity{ind(ii)} = 'Uncertain';
                    %                     end
                    
                    %% Delete events corresponding to the last channel
                    MarkingsValidated = HFOAnalysisResultsSingleInterval_ToPlot.MarkingsFR1;
                    MarkingsValidated = [MarkingsValidated,nEventValidity];
                    MarkingsValidated = [MarkingsValidated;Markings_ToAdd];
                    %                     indToDelete = ismember(MarkingsValidated(:,1),ListOfChannelsToPlot(end));
                    %                     MarkingsValidated(indToDelete,:) = [];
                    %                     nEventValidity(indToDelete,:) = [];
                    %                     strEventValidity(indToDelete,:) = [];
                    
                    %% Save variables
                    strChannelName = ElectrodeLabels{nChan};
                    strVariableName = sprintf(strFormat.SaveValidationResults,strEventType2,strrep(strrep(strFileName,'HFO_',''),'.mat',''),nGroup,nChan,strChannelName);
                    strVariablePath = [strValidationResultsFolder,strVariableName];
%                     save(strVariablePath,'MarkingsValidated','nChan','nEventValidity','strEventValidity','ListOfChannelsToPlot','strEventType','strChannelName','nGroup','nChannelRange_ToDelete','nChannelRange_ToDelete2')
                    
                else
                    %% Save variables
                    MarkingsValidated = [];
                    nEventValidity = [];
                    strEventValidity = {};
                    strChannelName = ElectrodeLabels{nChan};
                    strVariableName = sprintf(strFormat.SaveValidationResults,strEventType2,strrep(strrep(strFileName,'HFO_',''),'.mat',''),nGroup,nChan,strChannelName);
                    strVariablePath = [strValidationResultsFolder,strVariableName];
%                     save(strVariablePath,'MarkingsValidated','nChan','nEventValidity','strEventValidity','ListOfChannelsToPlot','strEventType','strChannelName','nGroup')
                end
            end
        else
            for nChan = 1:nNumberOfChannels
                %% Save variables
                MarkingsValidated = [];
                nEventValidity = [];
                strEventValidity = {};
                ListOfChannelsToPlot = ListOfElectrodesForDisplay(nChan,:);
                strChannelName = ElectrodeLabels{nChan};
                strVariableName = sprintf(strFormat.SaveValidationResults,strEventType2,strrep(strrep(strFileName,'HFO_',''),'.mat',''),nGroup,nChan,strChannelName);
                strVariablePath = [strValidationResultsFolder,strVariableName];
%                 save(strVariablePath,'MarkingsValidated','nChan','nEventValidity','strEventValidity','ListOfChannelsToPlot','strEventType','strChannelName','nGroup')
            end
        end
    end
end


