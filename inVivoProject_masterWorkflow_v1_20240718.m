% Master workflow to reproduce all analysis and plots included in Anja
% Payne's in vivo paper
% Written by Anja Payne
% Last modified: 07/29/2024

% Issues to resolve/currently working on:
%   - rerun step 5 when you have an hour to kill
%   - step 6: rate maps
%   - step 8: get the in-field spikes

% Steps:
%   1) Define pathway [done]
%   2) Read in excel files and create the main data structure [done]
%   3) Split the data into WT and KO [done]
%   4) Save the spike times in msec aligned to 0 [done]
%   5) Exclude spikes that occur during low velocity and format so that the
%      spikes are split by direction and trial number [done]
%   6) Get the rate maps [started, hold off]
%   7) Get spatial metrics and field barcode [started, hold off]
%   8) Get the in-field spikes
%   9) Get the theta modulation
%   10) Get the phase precession

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

% Settings: NA

% Inputs:
excelFolder = 'Z:\Anja\Data\In Vivo Data\AnimalInfo_ExcelFiles\'; 
excelFiles = {'AP_AnimalData_Cohort13_Track', 'AP_AnimalData_Cohort14_Track'}; 

% Outputs: 
mainDataStructure = getDataStructure_v1_20240718(excelFiles, excelFolder); 

%% Step 3: Split the file paths into WT and KO
clear;clc;

% Settings: 
fileNameBase = 'mainDataStructure'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '_settings.mat']);
mainDataStructureSettings = settings; % Empty file, for now

% Inputs: 
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
mainDataStructure = data;

% Outputs: 
filePaths = splitWTandKO_v1_20240722(mainDataStructure, mainDataStructureSettings); 

%% Step 4: Save the spike times (takes ~4 min)
clear;clc; tic; 

% Settings: 
fileNameBase = 'filePathsByGenotype';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '_settings.mat']);
filePathSettings = settings; % Empty file, for now

% Inputs: 
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
filePaths = data;

% Outputs: 
spikeTimes = getSpikeTimes_v1_20240725(filePaths, filePathSettings); toc

%% Step 5: Exclude spikes that occur during low velocity (takes ~50 min)
clear;clc;tic;

% Settings: 
fileNameBase = 'spikeTimes';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '_settings.mat']);
settings.velocity.threshold = 2; % velocity threshold = 2 cm/sec
settings.velocity.samplingRate = 16000; % sampling rate of video recording
settings.velocity.timeToAverage = 1; % time to average over in seconds
settings.rateMaps.trackWidth = 52; % cm; used to convert from pixels to cm
settings.rateMaps.trackLength = 80; % cm; used to convert from pixels to cm
settings.rateMaps.binSize = 4; % 4 cm bins

% Inputs: 
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
spikeTimes = data;

% Outputs: 
binnedSpikesByTrial = getHighVelocitySpikesByTrial_v1_20240725(spikeTimes, settings); toc

%% Step 6: Get the linearized rate maps
clear;clc;tic;

% Settings: 
fileNameBase = 'highVelocitySpikeTimes';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '_settings.mat']);
settings.rateMaps.trackSize = 264; 

% Input: 
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
binnedSpikesByTrial = data;

% Output: 
rateMaps = calculateRateMap_v1_20240718(binnedSpikesByTrial, settings); 

%% Step 7: Get spatial metrics and associated barcode
clear; clc; 

% Settings: 
fileNameBase = 'rateMaps';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '_settings.mat']);
rateMapSettings = settings; 
settings.rateMap.lowThresh = 0.1; 
settings.rateMap.highThresh = 0.5; 

% Inputs:
load([filePath{1}, '\', loadFileName, '.mat']);
rateMapStructure = data;

% Outputs:
data = getSpatialMetrics_v1_20240724(rateMapStructure, settings); 

%% Step 7: Get the in-field spikes for each field
clear;clc;tic;

% Settings: 
fileNameBase = 'highVelocitySpikeTimes';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '_settings.mat']);

% Inputs: 
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
binnedSpikes = data;

% Outputs: 

%% Step 9: Assign spikes and position to in-field
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
for iDir = 1:2; 
    %allSpikes{iAnimal}{iCluster}{iDir} = {}; 
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


