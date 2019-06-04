classdef Detections
    properties
        rParaPath
        frParaPath
        DataDir
        
        numDataFiles
        channelNames
        
        rHFOInfo
        frHFOInfo
        CoOccurenceInfo
    end
    
    methods
        function obj = setPDPaths(obj, rparaPath, frparaPath, dataDir)
            obj.rParaPath  = rparaPath;
            obj.frParaPath = frparaPath;
            obj.DataDir    = dataDir;
        end
        
        function obj = runDetector(obj, RefType, CondMulti, analDepth, channelContains, smoothBool, DetecType, ContThresh)
            if nargin < 2
                RefType = 'spec';
                CondMulti = false;
                analDepth = 3;
                channelContains = '';
                smoothBool = false;
                DetecType = 'RipAndFRip';
                ContThresh = 0.9;
            elseif nargin < 3
                CondMulti = false;
                analDepth = 3;
                channelContains = '';
                smoothBool = false;
                DetecType = 'RipAndFRip';
                ContThresh = 0.9;
            elseif nargin < 4
                analDepth = 3;
                channelContains = '';
                smoothBool = false;
                DetecType = 'RipAndFRip';
                ContThresh = 0.9;
            elseif nargin < 5
                channelContains = '';
                smoothBool = false;
                DetecType = 'RipAndFRip';
                ContThresh = 0.9;
            elseif nargin < 6
                smoothBool = false;
                DetecType = 'RipAndFRip';
                ContThresh = 0.9;
            elseif nargin < 7 
                DetecType = 'RipAndFRip';
                ContThresh = 0.9;
            elseif nargin < 8 
                ContThresh = 0.9;
            end
            % Paths
            rparaPath  = obj.rParaPath;
            frparaPath = obj.frParaPath;
            dataDir    = obj.DataDir;
            %% Ripples
            if isequal(DetecType,'Rip') || isequal(DetecType,'RipAndFRip')
                disp('Currently running Ripple detection.')
                obj.rHFOInfo = Detections.getClusterHFOData(...
                    rparaPath, dataDir, RefType, CondMulti, analDepth, channelContains, smoothBool);
            end  
            %% Fast Ripples
            if isequal(DetecType,'FRip') || isequal(DetecType,'RipAndFRip')
                disp('Currently running Fast Ripple detection.')
                obj.frHFOInfo = Detections.getClusterHFOData(...
                    frparaPath, dataDir, RefType,CondMulti, analDepth, channelContains, smoothBool);                        
            end
            %% Ripple And Fast Ripples
            if isequal(DetecType,'RipAndFRip')
                disp('Currently running Ripple and Fast Ripple overlap detection.')
                obj.CoOccurenceInfo = Detections.getClusterCoOccurenceInfo(obj.rHFOInfo, obj.frHFOInfo, ContThresh);%%%%%%%%%%%%%%%%%%
            end
            
            [~, obj.numDataFiles]  = Detections.getFileNames(dataDir, '*.mat');
            obj.channelNames   = obj.rHFOInfo.hfoInfo{1}.Data.channelNames; 
            
        end
        
        function [] = exportHFOsummary(obj)
            dataDir            = obj.DataDir;
            chanNames          = obj.channelNames;
            RipInfo            = obj.rHFOInfo;
            FRipInfo           = obj.frHFOInfo;
            CoOccurenceInf     = obj.CoOccurenceInfo;
            
            saveDir = [dataDir,'HFOSummary\'];
            if ~isdir(saveDir)
                mkdir(dataDir,'HFOSummary')
            end
            Detections.saveHFOCluster(RipInfo, FRipInfo, CoOccurenceInf, saveDir, chanNames)
        end
        
    end
    
    methods(Static)
        %% Processing HFO on mass   
        function clusterHFO = getClusterHFOData(paraPath, dataDir , RefType, CondMulti, analDepth, channelContains, smoothBool)
            [ListOfFiles, nbFiles] = Detections.getFileNames(dataDir, '*.mat');
            for iFile = 1:nbFiles
                disp(['Currently running interval: ',num2str(iFile), ' of ',num2str(nbFiles)])
                DataPath = [dataDir, ListOfFiles{iFile}];
                
                hfo = Detections.getHFOdata(paraPath, DataPath, RefType,...
                    CondMulti, analDepth, channelContains, smoothBool);
                if iFile == 1
                    nbChan   = length(hfo.baseline.maxNoisemuV);
                    NoiseMat = nan(nbFiles, nbChan);
                    BaseLMat = nan(nbFiles, nbChan);
                    HFOraMat = nan(nbFiles, nbChan);
                    
                    hfoInfo = cell(1,nbFiles);
                end
                NoiseMat(iFile, :) = hfo.baseline.maxNoisemuV;
                BaseLMat(iFile, :) = hfo.baseline.baselineThr;
                HFOraMat(iFile, :) = hfo.Events.Rates;
 
                % Strip Big variables of low information density
                hfo.Data = rmfield(hfo.Data,{'signal'});
                hfo.filtSig = [];
                hfo.baseline = rmfield(hfo.baseline,{'IndBaseline'});
                hfo.Events   = rmfield(hfo.Events,{'EventProp'});
                % and make it a struckt
                hfoInfo{iFile} = struct(hfo) ;
            end
            clusterHFO.channelNames = hfo.Data.channelNames;
            clusterHFO.numDataFiles = nbFiles;
            
            clusterHFO.Summary.NoiseThresholds = NoiseMat;
            clusterHFO.Summary.BaseLines       = BaseLMat;
            clusterHFO.Summary.HFOrates        = HFOraMat;
            
            clusterHFO.hfoInfo = hfoInfo;
        end
        % Processing a single HFO
        function hfo = getHFOdata(ParaPath, DataPath, RefType, CondMulti, analDepth, channelContains, smoothBool)
            if nargin < 3
                RefType = 'spec';
                CondMulti = false;
                analDepth = 3;
                channelContains = {''};
                smoothBool = false;
            elseif nargin < 4
                CondMulti = false;
                analDepth = 3;
                channelContains = {''};
                smoothBool = false;
            elseif nargin < 5
                analDepth = 3;
                channelContains = {''};
                smoothBool = false;
            elseif nargin < 6
                channelContains = {''};
                smoothBool = false;
            elseif nargin < 7
                smoothBool = false;
            end
            
            if  analDepth >= 1
                hfo = Core.HFO;
                hfo.ParaFileLocation = ParaPath;
                hfo.DataFileLocation = DataPath;
                hfo = getParaAndData(hfo, channelContains);
                hfo = getFilteredSignal(hfo, smoothBool); 
                if analDepth >= 2
                    hfo = getBasline(hfo);
                    if analDepth >= 3
                        hfo = getEvents(hfo, RefType, CondMulti);
                    end
                end
            end
        end
       
        %% co Occurence 
        function CoOccurence = getClusterCoOccurenceInfo(RhfoInfo, FRhfoInfo, ContThresh)
            nbFiles = RhfoInfo.numDataFiles;
            nbChan  = length(RhfoInfo.channelNames);
            
            EventInfo = cell(1,nbFiles);
            HFOrates = nan(nbFiles, nbChan);
            for iFile = 1:nbFiles
                rhfo = RhfoInfo.hfoInfo{iFile};
                frhfo = FRhfoInfo.hfoInfo{iFile};
%                 CoOc = Core.CoOccurence;
%                 EventInfo{iFile} = CoOc.runCoOccurence(rhfo, frhfo, ContThresh);
%                 EventInfo{iFile} = Detections.getECECoOccurence(rhfo, frhfo);
                EventInfo{iFile}   = Detections.getECECoOccurence(rhfo, frhfo);
                HFOrates(iFile, :) = EventInfo{iFile}.Rates.RippleANDFastRipple;
            end
            CoOccurence.EventInfo = EventInfo;
            CoOccurence.Summary.HFOrates = HFOrates;
        end

        %% Output
        % Export information
        function [] = saveHFOCluster(RipInfo, FRipInfo, CoOccurInfo, saveDir, chanNames)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rates %%%%%%%%%%%%%%%%%%%%
            RHFOraMat      = RipInfo.Summary.HFOrates;
            FRHFOraMat     = FRipInfo.Summary.HFOrates;
            RandFRHFOraMat = CoOccurInfo.Summary.HFOrates;
            [HFOArea, maskHFOArea, valMaxTHR ] = Detections.GetHFOAreaMat(RandFRHFOraMat');
             
            Detections.exportHFOrateDistribution(RHFOraMat, FRHFOraMat, RandFRHFOraMat, HFOArea, valMaxTHR, chanNames, saveDir)
            Detections.exportHFOrateDistribution([], [], RandFRHFOraMat, HFOArea, valMaxTHR, chanNames, saveDir,'R&FR')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Noise %%%%%%%%%%%%%%%%%%%%
            RNoiseMat   = RipInfo.Summary.NoiseThresholds;
            RBaseLMat   = RipInfo.Summary.BaseLines;
            
            FRNoiseMat  = FRipInfo.Summary.NoiseThresholds;
            FRBaseLMat  = FRipInfo.Summary.BaseLines;
            
            Detections.exportNoiseDist(RBaseLMat, FRBaseLMat, RNoiseMat, FRNoiseMat, HFOArea, chanNames, saveDir)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%
            % HFO intermediate inforamtion
            hfoSummary.Ripples     = RipInfo.hfoInfo;
            hfoSummary.FastRipples = FRipInfo.hfoInfo;
            hfoSummary.RandFR      = CoOccurInfo.EventInfo;
            % HFO summary inforamtion
            hfoSummary.Rates.RHFO        = RHFOraMat;
            hfoSummary.Rates.FRHFO       = FRHFOraMat;
            hfoSummary.Rates.RandFRHFO   = RandFRHFOraMat;
            
            hfoSummary.Area.HFOArea      = HFOArea;
            hfoSummary.Area.HFOthreshold = valMaxTHR;
            hfoSummary.Area.maskHFOArea  = maskHFOArea;
            
            hfoSummary.Noise.Ripples        = RNoiseMat;
            hfoSummary.Noise.FastRipples    = FRNoiseMat;
            hfoSummary.Baseline.Ripples     = RBaseLMat;
            hfoSummary.Baseline.FastRipples = FRBaseLMat;
            
            Detections.exportHFOdata(hfoSummary, saveDir)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reproducibiility
            
%             Detections.exportReproducibilityPlots(RandFRHFOraMat, HFOArea, maskHFOArea, chanNames, saveDir)
        end
              
        function [] = exportHFOdata(HFOSummaryMat, SaveDir)
            fileName = 'HFOSummaryMat';
            dataName = [SaveDir, fileName];
            save(dataName, fileName)
        end
        
        %% Visualizations
        % Noise and baseline
        function exportNoiseDist(RBaseLMat, FRBaseLMat, RNoiseMat,  FRNoiseMat, HFOArea, chanNames, saveDir)
            distPlot = Detections.plotDataDistBaseline(RBaseLMat, FRBaseLMat, RNoiseMat,  FRNoiseMat, HFOArea, chanNames, 'Baseline information');
            saveas(distPlot,[saveDir,'NoiseDistribution','.png'])
            close all
        end

        function fig = plotDataDistBaseline(RbaselineDataMat, FRbaselineDataMat, RnoiseDataMat, FRnoiseDataMat, HFOArea, chanNames, titleText)
            nbChannels = length(chanNames);
            fig = figure('units','normalized','outerposition',[0 0 1 1]);
            ax = Detections.makeFigureTight(gca, 6);
            
            
            hold on
            rbbar =bar(RbaselineDataMat',1,'blue');
            frbbar =bar(FRbaselineDataMat',0.75,'cyan');
            
            rnbar =bar(RnoiseDataMat',0.05,'green');
            frnbar =bar(FRnoiseDataMat',0.025,'yellow');
            hold off
            
            Ledg = legend([rbbar(1),frbbar(1),rnbar(1),frnbar(1)],{'ripple baseline','fast ripple baseline','ripple noise threshold','fast ripple noise threshold'});
            xtikkis = 1:nbChannels;
            
            set(gca,'xtick',[])
            ax2 = copyobj([gca,Ledg],gcf);
            ax2(1,1).XTickLabelRotation = 90;
            set(ax2(1,1),'ytick',[],'xtick',xtikkis(~HFOArea), 'xticklabel', chanNames(~HFOArea), 'fontsize', 8)
            ax3 = copyobj([gca,Ledg],gcf);
            ax3(1,1).XTickLabelRotation = 90;
            ax3(1,1).XTickLabelRotation = 90;
            set(ax3(1,1),'ytick',[],'xtick',xtikkis(HFOArea), 'XColor', 'red','xticklabel', chanNames(HFOArea), 'fontsize', 10)
            
            
            
            title(titleText,'fontsize', 20)
            annotation('textbox', [0.1, 0.86, 0.1, 0.1], 'string',  'muV')
        end
        
        % Event Rates
        function exportHFOrateDistribution(RHFOraMat, FRHFOraMat, RandFRHFOraMat, HFOArea, valMaxTHR, channelNames, saveDir, optionalStr)
            if nargin < 8
                optionalStr = '';
            end
            distPlot =  Detections.plotDataDistRates(RHFOraMat, FRHFOraMat, RandFRHFOraMat, HFOArea, valMaxTHR, channelNames,  'HFO rate distribution by channel');
            saveas(distPlot,[saveDir, optionalStr, 'HFORateDistribution','.png'])
            close all
        end         

        function fig = plotDataDistRates(RHFOraMat, FRHFOraMat, RandFRHFOraMat, HFOArea, valMaxTHR, chanNames, titleText)
            nbChannels = length(chanNames);
            fig = figure('units','normalized','outerposition',[0 0 1 1]);
            
            ylim([0 50])
            hold on
            if ~isempty(RHFOraMat) &&  isempty(FRHFOraMat) &&  isempty(RandFRHFOraMat) 
                rbar = bar(RHFOraMat', 0.3, 'blue');
                Ledg = legend(rbar(1),{'ripple'});
            end  
            if  isempty(RHFOraMat) && ~isempty(FRHFOraMat) &&  isempty(RandFRHFOraMat) 
                frbar = bar(FRHFOraMat', 0.6,  'cyan');
                Ledg = legend([frbar(1)],{'fast ripple'});
            end    
            if  isempty(RHFOraMat) &&  isempty(FRHFOraMat) && ~isempty(RandFRHFOraMat) 
                rfrbar = bar(RandFRHFOraMat', 2, 'red');
                Ledg = legend(rfrbar(1),{'ripple and fast ripple'});
                ylim([0 15])
            end
            if ~isempty(RHFOraMat) && ~isempty(FRHFOraMat) && ~isempty(RandFRHFOraMat)
                rbar   = bar(RHFOraMat', 0.3, 'blue');
                frbar  = bar(FRHFOraMat', 0.6,  'cyan');
                rfrbar = bar(RandFRHFOraMat', 2, 'red');
                Ledg   = legend([rbar(1), frbar(1), rfrbar(1)],{'ripple','fast ripple','ripple and fast ripple'});
            end
            hold off
            
            xtikkis = 1:nbChannels;
            
            threshLine = line([0, nbChannels],[valMaxTHR, valMaxTHR]);
            threshLine.LineWidth = 2;
            
            ax = Detections.makeFigureTight(gca, 6);
            
            %%%%%%%%%%%%%%%%
            set(gca,'xtick',[])
            ax2 = copyobj([gca,Ledg],gcf);
            ax2(1,1).XTickLabelRotation = 90;
            set(ax2(1,1),'ytick',[],'xtick',xtikkis(~HFOArea), 'xticklabel', chanNames(~HFOArea), 'fontsize', 8)
            ax3 = copyobj([gca,Ledg],gcf);
            ax3(1,1).XTickLabelRotation = 90;
            set(ax3(1,1),'ytick',[],'xtick',xtikkis(HFOArea), 'XColor', 'red','xticklabel', chanNames(HFOArea), 'fontsize', 10)
            %%%%%%%%%%%%%%%%%%%%%

            
            
            title(titleText,'fontsize', 20)
            annotation('textbox', [0.1, 0.86, 0.1, 0.1], 'string',  'Event/min')
        end
        
        %Reproducibility
        function exportReproducibilityPlots(RandFRHFOraMat, HFOArea, maskHFOArea, chanNames, saveDir)
            
            [RandomDistofSP, ActualDistofSP] = Misc.Reproducibility.getTestReTest(RandFRHFOraMat, 10000);
            
            close
            custom_map = [1,1,1; 0.8,0.8,0.8; 1,0,0];
            
            distScalarFig = Detections.plotdistScalarFig(RandomDistofSP, ActualDistofSP, custom_map);
            saveas(distScalarFig, [saveDir, 'ReproducibilityEstimate', '.png'])
            close

            figmaskHFO = Detections.plotFigMaskHFO(chanNames, maskHFOArea, HFOArea, custom_map);
            saveas( figmaskHFO, [saveDir, 'HFOArea', '.png'])
            close
        end
        
        function distScalarFig = plotdistScalarFig(RandomDistofSP, ActualDistofSP, custom_map)
            PercentileVal = prctile(RandomDistofSP(:), 95);
            colormap(custom_map)
            binStarts = 0:0.01:1;
            distScalarFig = figure('units','normalized','outerposition',[0 0 1 1]);
            hold on
            hist1 = histogram(RandomDistofSP(:), binStarts,'Normalization','probability','facecolor',[0.8,0.8,0.8],'facealpha', 0.85);
            hist2 = histogram(ActualDistofSP(:), binStarts,'Normalization','probability','facecolor',[1,0,0],'facealpha', 0.85);
            PercentileLine = line([PercentileVal, PercentileVal], [0, 0.3]);
            PercentileLine.LineWidth = 2;
            hold off
            ylabel('Probability Desnity')
            xlabel('Scalar Product')
            title('Random HFO rate distribution in time vs. Measured HFO rate distribution in time.')
            legend([hist1(1), hist2(1), PercentileLine(1)],{'Random Distribution','Measured Distibution','95th percentile of random distribution.'});
            ax = Detections.makeFigureTight(gca, 2);
        end
        
        function figmaskHFO = plotFigMaskHFO(chanNames,maskHFOArea, HFOArea, custom_map)
            nbChannels = length(chanNames);
            figmaskHFO = figure('units','normalized','outerposition',[0 0 1 1]);
            
            ax = Detections.makeFigureTight(gca, 3);
            
            nbIntevals = size(maskHFOArea,2);
            TempSeq = num2cell(1:nbIntevals);
            TempStrSeq = cellfun(@num2str, TempSeq,'un',0);
            xtics = [TempStrSeq,'Blank','Average'];
            
            imagesc([maskHFOArea, -ones(length(HFOArea),1), HFOArea])
            xticklabels(xtics)
            xticks(1:(nbIntevals+2))
             
            colormap(custom_map)
            ytikkis = 1:nbChannels;
            set(gca,'ytick',[])
            ax2 = copyobj(gca,gcf);
            set(ax2(1,1), 'xtick', [], 'ytick', ytikkis(~HFOArea), 'yticklabel', chanNames(~HFOArea), 'fontsize', 6)
            ax3 = copyobj(gca,gcf);
            

            set(ax3(1,1), 'xtick', [], 'xticklabel', xtics, 'ytick', ytikkis(HFOArea), 'YColor', 'red','yticklabel', chanNames(HFOArea), 'fontsize', 9)
            
            xlabel('Interval')
            ylabel('Channel')
            title('Channel inclusion in HFO area.')
        end
        
        function ax = makeFigureTight(gca, tightness)
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset*tightness;
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
        end
        
        
        %% HFO area
        function [HFOArea, maskHFOArea,valMaxTHRval ] = GetHFOAreaMat(RandFRHFOraMat)
            % Within interval 95th pecentile selection for HFO area(channel space)
            % Hence, this code answers which channels in a given interval
            % of recoding contribute the most HFO
            valMaxTHR = prctile(RandFRHFOraMat,95);
            valMaxTHRval = max(valMaxTHR);
            valMaxTHR = repmat(valMaxTHR, size(RandFRHFOraMat, 1), 1);
            maskNoHFOInterval = (valMaxTHR == 0);
            maskHFOArea =  (RandFRHFOraMat >= valMaxTHR);
            maskHFOArea(maskNoHFOInterval) = 0;
            
            % Aggregation accross intervals
            RandFRHFOraVec = mean(maskHFOArea ,2);
            % accros aggregate  95th percentile section for HFO area(channel space)
            % Hence, this code answers which channels across all intervals
            % of recoding contribute the most consistently HFO
            [HFOArea, ~] = Detections.GetHFOAreaVec(RandFRHFOraVec);
        end
        
        function [ HFOArea, valMaxTHRval ] = GetHFOAreaVec(RandFRHFOraVec)
            % Within interval 95th pecentile selection for HFO area(channel space)
            % Hence, this code answers which channels in a given interval
            % of recoding contribute the most HFO
            valMaxTHR = prctile(RandFRHFOraVec,95);
            valMaxTHRval = max(valMaxTHR);
            valMaxTHR = repmat(valMaxTHR, size(RandFRHFOraVec, 1), 1);
            HFOArea =  (RandFRHFOraVec >= valMaxTHR);
        end
        
        %% File Utility        
        function [ListOfFiles, nbFiles] = getFileNames(LoadDirectory, extension)
            % This function returns the names of .mat files in a given directory.
            % Input: directory location as string e.g. '/home/andries/DataForProjects/SleepData/Patient1/'
            % Output: cell with string entries e.g.  {'pat_1_night_1_interval_1.mat'}...                                             .
            
            if 7~=exist(LoadDirectory,'dir')
                error('This directory does not exist.')
            end
            
            addpath(genpath(LoadDirectory))
            ListOfFiles = dir([LoadDirectory, extension]);
            ListOfFiles = {ListOfFiles.name}';
            nbFiles = length(ListOfFiles);
        end

        %% THis is the old stuff you have to use to mathc up with the old data
        function CoOccurenceInfo  = getECECoOccurence(Rhfo, FRhfo)
            nbChan = Rhfo.Data.nbChannels;
            CoOccurenceInfo.maskCell.Ripples             = cell(nbChan,1);
            CoOccurenceInfo.maskCell.FastRipples         = cell(nbChan,1);
            CoOccurenceInfo.maskCell.RippleANDFastRipple = cell(nbChan,1);
            merged = cell(nbChan,1);
            for iChan = 1:nbChan
                
                maskOfContainmentRFR = Detections.getIndOfContainment(Rhfo, FRhfo, iChan); %Ind (mark = 3)
                % ~maskOfContainmentRFR % is (mark = 1)        %         (mark1 | mark3)  (OLD RIPPLE MASK)
                maskOfContainmentFRR = Detections.getIndOfContainment(FRhfo, Rhfo, iChan); %Ind (OLD FAST RIPPLE MASK)
                
                CoOccurenceInfo.maskCell.Ripples{iChan}             =  true(1,length(maskOfContainmentRFR));
                CoOccurenceInfo.maskCell.FastRipples{iChan}         =  {maskOfContainmentRFR, ~maskOfContainmentFRR}; % for a concatinated Ripples and fast ripples
                CoOccurenceInfo.maskCell.RippleANDFastRipple{iChan} =  maskOfContainmentRFR; % taken from the ripples
                
                merged{iChan} = [maskOfContainmentRFR, ~maskOfContainmentFRR];
                
            end
            CoOccurenceInfo.Count.Ripple = cellfun(@sum,CoOccurenceInfo.maskCell.Ripples);
            CoOccurenceInfo.Count.FastRipple = cellfun(@sum,merged);
            CoOccurenceInfo.Count.RippleANDFastRipple = cellfun(@sum,CoOccurenceInfo.maskCell.RippleANDFastRipple);
            
            durationMin = (Rhfo.Data.nbSamples/Rhfo.Data.sampFreq)/60;
            
            CoOccurenceInfo.Rates.Ripple              = CoOccurenceInfo.Count.Ripple/durationMin;
            CoOccurenceInfo.Rates.FastRipple          = CoOccurenceInfo.Count.FastRipple/durationMin;
            CoOccurenceInfo.Rates.RippleANDFastRipple = CoOccurenceInfo.Count.RippleANDFastRipple/durationMin;
        end

        function IndOfContainment = getIndOfContainment(Rhfo, FRhfo, iChan)
            fs      = Rhfo.Data.sampFreq;
            nbEvent = Rhfo.Events.EventNumber(iChan);
            Ind1    = zeros(1,nbEvent);
            for iEvent = 1:nbEvent
                Rs  = Rhfo.Events.Markings.start{iChan}/fs;
                Re  = Rhfo.Events.Markings.end{iChan}/fs;
                FRe = FRhfo.Events.Markings.end{iChan}/fs;
                FRs = FRhfo.Events.Markings.start{iChan}/fs;
                
                maskRipFRipStart = (Rs(iEvent) < FRe);
                maskRipFRipEnd   = (Re(iEvent) > FRs);
                maskRipFRipBoth  = (maskRipFRipStart & maskRipFRipEnd);
                
                Ind1(iEvent)     =  any(maskRipFRipBoth);
            end
            IndOfContainment = logical(Ind1);
        end
        
    end
end
   