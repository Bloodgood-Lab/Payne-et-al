% Master workflow to reproduce all analysis and plots included in Anja
% Payne's in vivo paper
% Written by Anja Payne
% Last modified: 07/29/2024

% Issues to resolve/currently working on:
%   - I think the settings should be in the main data structure instead of
%     a separate file
%   - i feel like the processed data should be in folders... maybe each
%     folder name should have the last save date appended?
%   - step 9: split into high and low firing rate cells

% Steps:
%   1) Define pathway [done]
%   2) Read in excel files and create the main data structure [done]
%   3) Split the data into WT and KO [done]
%   4) Save the spike times in msec aligned to 0 [done]
%   5) Exclude spikes that occur during low velocity and format so that the
%      spikes are split by direction and trial number [done]
%   6) Get the rate maps [done]
%   7) Using firing rates, divide data into high-firing (putative place
%      cells) and low-firing [done]
%   8) Get spatial metrics and field barcode [done]

%   9) Get the in-field spikes
%   10) Additional rate map analysis on the high-firing and low-firing and
%       in-field firing
%   11) Placeholder for all the stability analysis
%   12) Get the theta modulation
%   13) Get the phase precession

% Notes about this code:
%   - When you choose to save processed data (as prompted by the pop-up
%     window(s), the .mat file will always be saved with increasing 
%     versions and with the data appended.
%   - Each step will give you the option of saving the output. The next 
%     step will look for that saved file so it's a good idea to save as you
%     go. The code will save two files, one with the data and one with the
%     settings. 

%% Step 1: Define pathways
addpath(genpath('Z:\Anja\Paper\Matlab Code')); 

%% Step 2: Read in excel files and create the main data structure
clear;clc;

% Inputs:
excelFolder = 'Z:\Anja\Data\In Vivo Data\AnimalInfo_ExcelFiles\'; 
excelFiles = {'AP_AnimalData_Cohort13_Track', 'AP_AnimalData_Cohort14_Track'}; 

% Outputs: 
mainDataStructure = getDataStructure_v1_20240718(excelFiles, excelFolder); 

%% Step 3: Split the file paths into WT and KO
clear;clc;

% Inputs: 
fileNameBase = 'mainDataStructure'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
mainDataStructure = data;

% Settings
mainDataStructureSettings = struct(); 

% Outputs: 
filePaths = splitWTandKO_v1_20240722(mainDataStructure, mainDataStructureSettings); 

%% Step 4: Save the spike times (takes ~4 min)
clear;clc; tic; 

% Inputs: 
fileNameBase = 'filePathsByGenotype';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
filePaths = data; 

% Settings: 
filePathSettings = settings;

% Outputs: 
spikeTimes = getSpikeTimes_v1_20240725(filePaths, filePathSettings); toc

%% Step 5: Exclude spikes that occur during low velocity (takes ~55 min)
clear;clc;tic;

% Inputs: 
fileNameBase = 'spikeTimes';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
spikeTimes = data; spikeTimeSettings = settings; 

% Settings: 
spikeTimeSettings = settings;
spikeTimeSettings.velocity.threshold = 2; % velocity threshold = 2 cm/sec
spikeTimeSettings.velocity.samplingRate = 16000; % sampling rate of video recording
spikeTimeSettings.velocity.timeToAverage = 1; % time to average over in seconds
spikeTimeSettings.rateMaps.trackWidth = 52; % cm; used to convert from pixels to cm
spikeTimeSettings.rateMaps.trackLength = 80; % cm; used to convert from pixels to cm
spikeTimeSettings.rateMaps.binSize = 4; % 4 cm bins
[suppDataPath, ~, ~] = fileparts(filePath{1}); % location to store binned position files

% Outputs: 
binnedSpikesByTrial = getHighVelocitySpikesByTrial_v1_20240725(spikeTimes, spikeTimeSettings, suppDataPath); toc

%% Step 6: Get the linearized rate maps (takes ~3 min)
clear;clc;tic;

% Input: 
fileNameBase = 'highVelocitySpikeTimes';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
binnedSpikesByTrial = data; binnedSpikesByTrialSettings = settings; 

% Settings: 
binnedSpikesByTrialSettings.rateMaps.trackSize = 264; 

% Output: 
rateMaps = calculateRateMap_v1_20240718(binnedSpikesByTrial, binnedSpikesByTrialSettings); toc

%% Step 7: Split into high-firing and low-firing cells (takes seconds)
clear;clc;tic;

% Inputs:
fileNameBase = 'rateMaps';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '.mat']);
rateMaps = data; rateMapSettings = settings; 

% Settings: 
rateMapSettings.firingRates.meanThresh = 0.1; 
rateMapSettings.firingRates.maxThresh = 1; 

% Outputs:
outputData = splitLowAndHighFR_v01_20240802(rateMaps, rateMapSettings); toc

%% Step 8: Get spatial metrics and associated barcode (takes seconds)
clear;clc;tic;

% Inputs:
fileNameBase = 'rateMapsByFiringRate';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '.mat']);
rateMapByFRStructure = data; rateMapByFRSettings = settings; 

% Settings: 
rateMapByFRSettings.rateMaps.lowThresh = 0.1;  % fields will be counted as contiguous bins above 10% of max
rateMapByFRSettings.rateMaps.highThresh = 0.5; % fields must include one bin that is above 50% of the max

% Outputs:
spatialMetrics = getSpatialMetrics_v1_20240724(rateMapByFRStructure, rateMapByFRSettings); toc
    
%% Step 9: Get the in-field spikes for each field
clear;clc;tic;

% Inputs: 
fileNameBase = 'spatialMetrics';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
binnedSpikes = data; binnedSpikeSettings = settings;

% Settings: 

% Outputs: 

%% Step: Assign spikes and position to in-field
% Pre-allocate based on the length of the data
% if length(cw_spike_bins) > length(ccw_spike_bins); 
%     spike_logical{iAnimal}{iCluster} = zeros(length(cw_spike_bins),2);
%     allSpikes{iAnimal}{iCluster} = NaN(length(ccw_spike_bins),2);             
% elseif length(ccw_spike_bins) > length(cw_spike_bins); 
%     spike_logical{iAnimal}{iCluster} = zeros(length(ccw_spike_bins),2);
%     allSpikes{iAnimal}{iCluster} = NaN(length(ccw_spike_bins),2);             
% end
% if length(cw_pos) > length(ccw_pos); 
%     pos_logical{iAnimal}{iCluster} = zeros(length(cw_pos),2);
% elseif length(ccw_pos) > length(cw_pos); 
%     pos_logical{iAnimal}{iCluster} = zeros(length(ccw_pos),2);
% end
    for iGenotype = 1:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype}); 
        for iAnimal = 1:length(genotypeData); 
            if isempty(genotypeData{iAnimal}) == 1; 
                continue
            else
                [~,n] = size(genotypeData{iAnimal});
                for iCluster = 1:n;
                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                    for iDir = 1:2; 
                        if iDir == 1; 
                            %allSpikes{iAnimal}{iCluster}{iDir} = {}; 
                            barcode = 
                            for iField = 1:max(barcode{iAnimal}{iCluster, iDir}); 
                                iBinField = find(barcode{iAnimal}{iCluster, iDir} == iField); 
                                if iDir == 1; 
                                    % Determine which positions occur in the 'in-field' bin
                                    binIndex_pos = ismember(cw_pos, iBinField); 
                                    % Determine which spikes occur in the 'in-field' bin
                                    binIndex_spikes = ismember(cw_spike_bins, iBinField); 
                                    allSpikes{iAnimal}{iCluster}(1:length(cw_spikes),iDir) = cw_spikes;                            
                                elseif iDir == 2; 
                                    % Determine which positions occur in the 'in-field' bin                    
                                    binIndex_pos = ismember(ccw_pos, iBinField); 
                                    % Determine which spikes occur in the 'in-field' bin
                                    binIndex_spikes = ismember(ccw_spike_bins, iBinField);  
                                    allSpikes{iAnimal}{iCluster}(1:length(ccw_spikes),iDir) = ccw_spikes;                                                             
                                end
                                [m,~] = size(allSpikes{iAnimal}{iCluster});
                                pos_logical{iAnimal}{iCluster}(binIndex_pos, iDir) = iField; 
                                spike_logical{iAnimal}{iCluster}(1:m, iDir) = 0;  
                                spike_logical{iAnimal}{iCluster}(binIndex_spikes, iDir) = iField;  
                            end
                        end
                    end
                end
            end
        end
    end

%% Step 8: (Placeholder for stability stuff)

