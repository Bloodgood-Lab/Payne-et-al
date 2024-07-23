% Master workflow to reproduce all analysis and plots included in Anja
% Payne's in vivo paper
% Written by Anja Payne
% Last modified: 07/18/2024

% Issues to resolve:
%   - For step 3, set it up so that the input is always the last saved
%     version of the mainDataStructure
%   - Implement step 4: get linearized rate maps
%   - Clean up step 3, right now it's sort of acting as a placeholder for
%     what the output would be so that that can be fed into step 4
%   - Implement step 1
%   - Implement step 2

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

%% Step 3: Separate into low-firing and high-firing cells
% Settings:
maxFRcutOff     = 1; 
meanFRcutOff    = 0.1;
% Input:

% Output:

%% Step 4: Get the linearized rate maps
% Settings: 
generatePlots = 0; % set to 0 if you don't want the plots and 1 if you do
% Input: 
% Output: 







