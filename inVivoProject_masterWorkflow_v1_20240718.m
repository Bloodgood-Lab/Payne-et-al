% Master workflow to reproduce all analysis and plots included in Anja
% Payne's in vivo paper
% Written by Anja Payne
% Last modified: 07/18/2024

% Issues to resolve/currently working on:
%   - spike times: settings does not yet exist. At what point should I
%     create it?

% Steps:
%   1) Define pathway [done]
%   2) Read in excel files and create the main data structure [done]
%   3) Split the data into WT and KO [done]
%   4) Save the spike times in msec aligned to 0
%   4) Exclude spikes that occur during low velocity [started, hold off]
%   5) Get the rate maps [started, hold off]
%   6) Get spatial metrics and field barcode [current WIP]
%   7) Get the in-field spikes
%   8) Get the theta modulation
%   9) Get the phase precession

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

%% Step 4: Save the spike times (takes ~3 min)
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



%% Step 4: 
clear;clc;

% Settings: 
fileNameBase = 'filePathsByGenotype';
trackWidth = 52; trackLength = 80; % In cm
track{1} = trackWidth; track{2} = trackLength; 

% Inputs: 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
filePaths = data;

% Outputs: 
directoryPath = ''; 
for iGenotype = 1%:length(fieldnames(data));
    genotypes = fieldnames(data); 
    genotypeData = data.(genotypes{iGenotype}); 
    for iAnimal = 1%:length(genotypeData); 
        for iCluster = 1:2%:length(genotypeData{iAnimal});
            display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
            %% Step 1: Load the data
            % If this data has been loaded in a previous iteration,
            % skip loading it now
            newDirectoryPath = genotypeData{iAnimal}(iCluster).directory;
            if strcmp(directoryPath, newDirectoryPath) == 1;
                % Do nothing
            elseif strcmp(directoryPath, newDirectoryPath) == 0;
                directoryPath = newDirectoryPath; 
                % Load the position information
                idcs = strfind(directoryPath{1}, '\');
                parent_directory = directoryPath{1}(1:idcs(end)-1);
                load([parent_directory, '\Position.mat']);
                % Load the CSC data
                CSC_file = [directoryPath{1}, '\CSC1.ncs']; 
                csc_data = readCSC(CSC_file); 
                ephys_t0 = csc_data.ts(1); 
            end
            % Load the spike data for the relevant cluster
            spikeDataPath = strcat(directoryPath, '\', genotypeData{iAnimal}(iCluster).fileName, '.mat');
            load(spikeDataPath{1});
            spike_time_msec = TS/10;
            
            x_pos = aligned_pos(:,1); 
            y_pos = aligned_pos(:,2);
            widthConversionFactor = track{1}/(max(x_pos) - min(x_pos)); 
            lengthConversionFactor = track{2}/(max(y_pos) - min(y_pos)); 
            x_pos_cm = x_pos.*widthConversionFactor;
            y_pos_cm = y_pos.*lengthConversionFactor;
                
            %% Step 2: Exclude spikes during low velocity
            [include_spikes, ave_vel, std_vel, vel_thresh, velocity] = excludeLowVelocity(x_pos_cm, y_pos_cm, t_clipped, spike_time_msec, ephys_t0); 

                
        end
    end
end

%% Step 5: Get the linearized rate maps
clear; clc;

% Settings: 
fileNameBase = 'filePathsByGenotype';
trackWidth = 52; trackLength = 80; % In cm
track{1} = trackWidth; track{2} = trackLength; 
%generatePlots = 0; % set to 0 if you don't want the plots and 1 if you do

% Input: 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
filePaths = data;

% Output: 
clc;
calculateRateMap_v1_20240718(filePaths, track); 

%% Step 6: Get spatial metrics and field barcode
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
clear; clc; 

% Settings: 
fileNameBase = 'rateMaps';
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '_settings.mat']);
rateMapSettings = settings; 
settings.rateMap.lowThresh = 0.1; 
settings.rateMap.highThresh = 0.5; 

