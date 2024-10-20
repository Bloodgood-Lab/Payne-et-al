% Master workflow to reproduce all analysis and plots included in Anja
% Payne's in vivo paper
% Written by Anja Payne
% Last modified: 08/15/2024

% Steps:
%   1) Define pathway
%   2) Read in excel files and create the main data structure 
%   3) Split the data into WT and KO
%   4) Save the spike times in msec aligned to 0 
%   5) Exclude spikes that occur during low velocity and format so that the
%      spikes are split by direction and trial number
%   6) Get the rate maps
%   7) Using firing rates, divide data into high-firing (putative place
%      cells) and low-firing
%   8) Get spatial metrics and field barcode
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
[processedDataFolder, ~, ~] = fileparts(filePath{1});
mainDataStructure = data;

% Settings
mainDataStructureSettings = struct(); 

% Outputs: 
filePaths = splitWTandKO_v1_20240722(mainDataStructure, mainDataStructureSettings, processedDataFolder); 

%% Step 4: Save the spike times (takes ~4 min)
clear;clc; tic; 

% Inputs: 
fileNameBase = 'filePathsByGenotype';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
filePaths = data; 

% Settings: 
filePathSettings = settings;

% Outputs: 
spikeTimes = getSpikeTimes_v1_20240725(filePaths, filePathSettings, processedDataFolder); toc

%% Step 5: Exclude spikes that occur during low velocity (takes ~55 min)
clear;clc;tic;

% Inputs: 
fileNameBase = 'spikeTimes';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
spikeTimes = data; spikeTimeSettings = settings; 

% Settings: 
spikeTimeSettings = settings;
spikeTimeSettings.velocity.threshold = 2; % velocity threshold = 2 cm/sec
spikeTimeSettings.velocity.samplingRate = 16000; % sampling rate of video recording
spikeTimeSettings.velocity.timeToAverage = 1; % time to average over in seconds when calculating velocity
spikeTimeSettings.rateMaps.trackWidth = 52; % cm; used to convert from pixels to cm
spikeTimeSettings.rateMaps.trackLength = 80; % cm; used to convert from pixels to cm
spikeTimeSettings.rateMaps.binSize = 4; % 4 cm bins

% Outputs: 
binnedSpikesByTrial = getHighVelocitySpikesByTrial_v1_20240725(spikeTimes, spikeTimeSettings, processedDataFolder); toc

%% Step 6: Get the linearized rate maps (takes ~3 min)
clear;clc;tic;

% Input: 
fileNameBase = 'highVelocitySpikeTimes';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
binnedSpikesByTrial = data; binnedSpikesByTrialSettings = settings; 

% Settings: 
binnedSpikesByTrialSettings.rateMaps.trackSize = 264; 

% Output: 
rateMaps = calculateRateMap_v1_20240718(binnedSpikesByTrial, binnedSpikesByTrialSettings, processedDataFolder); toc

%% Step 7: Split into high-firing and low-firing cells (takes seconds)
clear;clc;tic;

% Inputs:
fileNameBase = 'rateMaps';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '.mat']);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
rateMaps = data; rateMapSettings = settings; 

% Settings: 
rateMapSettings.firingRates.meanThresh = 0.1; 
rateMapSettings.firingRates.maxThresh = 1; 

% Outputs:
dataByFiringRate = splitLowAndHighFR_v01_20240802(rateMaps, rateMapSettings, processedDataFolder); toc

%% Step 8: Get spatial metrics and associated barcode (takes seconds)
clear;clc;tic;

% Inputs:
fileNameBase = 'rateMapsByFiringRate';
folderMessage = 'Select directory with data to analyze'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase, folderMessage);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '.mat']);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
rateMapByFRStructure = data; rateMapByFRSettings = settings; 

% Settings: 
rateMapByFRSettings.rateMaps.lowThresh = 0.1;  % fields will be counted as contiguous bins above 10% of max
rateMapByFRSettings.rateMaps.highThresh = 0.5; % fields must include one bin that is above 50% of the max
rateMapByFRSettings.rateMaps.biggerFieldModification = 2; % for use in phase precession control analysis: extend field by X bins
rateMapByFRSettings.rateMaps.smallerFieldModification = -2; % for use in phase precession control analysis: reduce field by X bins

% Outputs:
spatialMetrics = getSpatialMetrics_v1_20240724(rateMapByFRStructure, rateMapByFRSettings, processedDataFolder); toc
    
%% Step 9: Get the in-field spikes for each field (takes seconds)
clear;clc;tic;

% Inputs: 
fileNameBase = 'spatialMetrics';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
binnedSpikes = data; binnedSpikeSettings = settings;

% Settings: 

% Outputs: 
inFieldSpkTimes = getInFieldSpikes_v1_20240805(binnedSpikes, binnedSpikeSettings, processedDataFolder); toc

%% Step 10: Additional rate map analysis

%% Step 11: Stability analysis

%% Step 12: Get the theta modulation (takes ~1.5 hours)
clear;clc;tic;
display(['Estimated time to finish is ', datestr(datetime('now')+hours(1.5), 'HH:MM:SS')]); 

% Inputs: 
fileNameBase = 'inFieldSpkTimes';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
inFieldSpkTimes = data; inFieldSpkTimesSettings = settings;

% Settings: 
inFieldSpkTimesSettings.theta.frequencyBand = [4,12]; 

% Outputs: 
thetaData = getThetaModulation_v1_20240806(inFieldSpkTimes, inFieldSpkTimesSettings, processedDataFolder); toc

%% Step 13: Get the phase precession (takes ~8 min)
clear;clc;close all; tic;

% Inputs: 
fileNameBase = 'theta';
folderMessage = 'Select directory with data to analyze'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase, folderMessage);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
thetaData = data; thetaSettings = settings;

% Settings: 
thetaSettings.phasePrecession.spatialBinThreshold = 0; % minimum number of spatial bins needed
thetaSettings.phasePrecession.spikeThreshold = 5; % minimum number of spikes needed
thetaSettings.phasePrecession.slopeRange = [-4:0.3:4]; % range of slopes to try to fit 
thetaSettings.phasePrecession.significanceThreshold = 1; % maximum acceptable significance level of line fit
thetaSettings.phasePrecession.trialThreshold = 5; % minimum number of trials
thetaSettings.phasePrecession.ISIthreshold = 1000; % max time between spikes in msec
thetaSettings.phasePrecession.fieldsToAnalyze = 'all fields'; % Could also be 'best field'
thetaSettings.phasePrecession.positionType = 'binned';
thetaSettings.phasePrecession.plot = 'yes'; % Determines which plots will be generated; 'yes' plots all while 'relationshipsOnly' only plots population data
thetaSettings.phasePrecession.normalized = 'yes'; % Is the field normalized?
thetaSettings.phasePrecession.circularity = 'none'; % How is circularity accounted for? Could also be set to 'shift'
thetaSettings.phasePrecession.fit = 'circularSlope'; % Using a linear or circular fit?
thetaSettings.phasePrecession.timeRange = 5*125; % Over what range of time should spikes occur? 625 msec = 5 theta cycles

% Outputs: 
[phasePrecessionData, phasePrecessionSettings] = getPhasePrecession_v1_20240806(thetaData, thetaSettings, processedDataFolder); toc 
%%
clc;
phasePrecessionSettings.phasePrecession.plot = 'shuffles'; 
plotPhasePrecession_v1_20240827(phasePrecessionData, phasePrecessionSettings); 

%% Step 13B: Control analyses related to the phase precession
clear;clc;close all; tic;

% Inputs: 
fileNameBase = 'phasePrecession';
folderMessage = 'Select directory with data to analyze'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase, folderMessage);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
phasePrecessionData = data; phasePrecessionSettings = settings;


% Settings: 


% Outputs: 
[data, settings] = getPhasePrecessionWithChangingFieldSize_v1_20240924(phasePrecessionData, phasePrecessionSettings, processedDataFolder)

%%





