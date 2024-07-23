% Master workflow to reproduce all analysis and plots included in Anja
% Payne's in vivo paper
% Written by Anja Payne
% Last modified: 07/18/2024

% Issues to resolve/currently working on:
%   - Save the output from step 3 (split into WT and KO)

% Steps:
%   1) Define pathway
%   2) Read in excel files and create the main data structure
%   3) Split the data into WT and KO
%   4) Get spiking metrics (i.e. firing rate and burst index)
%   5) Use the spiking metrics to separate into low firing cells and high
%      firing cells
%   6) Get the linearized rate maps

% Notes about this code:
%   - When you choose to save processed data (as prompted by the pop-up
%     window(s), the .mat file will always be saved with increasing 
%     versions and with the data appended.

%% Step 1: Define pathways
addpath(genpath('Z:\Anja\Paper\Matlab Code')); 

%% Step 2: Read in excel files and create the main data structure
% Settings: NA

% Inputs:
excelFolder = 'Z:\Anja\Data\In Vivo Data\AnimalInfo_ExcelFiles\'; 
excelFiles = {'AP_AnimalData_Cohort13_Track', 'AP_AnimalData_Cohort14_Track'}; 

% Outputs: 
mainDataStructure = getDataStructure_v1_20240718(excelFiles, excelFolder); 

%% Step 3: Split the data into WT and KO
% Settings: 
fileNameBase = 'mainDataStructure'; 

% Inputs: 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
mainDataStructure = data;

% Outputs: 
filePaths = splitWTandKO_v1_20240722(mainDataStructure); 

%% Step 4: Get the linearized rate maps
% Settings: 
generatePlots = 0; % set to 0 if you don't want the plots and 1 if you do
% Input: 
% Output: 







