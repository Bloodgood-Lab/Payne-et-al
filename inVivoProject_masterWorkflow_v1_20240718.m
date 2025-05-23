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
highVelocitySpikeTimes = getHighVelocitySpikesByTrial_v1_20240725(spikeTimes, spikeTimeSettings, processedDataFolder); toc

%% Step 6: Look at spike timing (bursts vs. singles)
clear;clc;tic;

% Inputs: 
fileNameBase = 'highVelocitySpikeTimes';
folderMessage = 'Select directory with data to analyze'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase, folderMessage);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
highVelSpkTimes = data; highVelSpkTimesSettings = settings; 

% Settings: 

% Outputs: 
spikeTiming = getSpikeTimingData_v1_20241124(highVelSpkTimes, highVelSpkTimesSettings, processedDataFolder); toc

%% Step 7: Get the linearized rate maps (takes ~3 min)
clear;clc;tic;

% Input: 
fileNameBase = 'spikeTiming';
folderMessage = 'Select directory with data to analyze'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase, folderMessage);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
spikeTimingData = data; spikeTimingSettings = settings; 

% Settings: 
spikeTimingSettings.rateMaps.trackSize = 264; 

% Output: 
rateMaps = calculateRateMap_v1_20240718(spikeTimingData, spikeTimingSettings, processedDataFolder); toc

%% Step 8: Split into high-firing and low-firing cells (takes seconds)
clear;clc;tic;

% Inputs:
fileNameBase = 'rateMaps';
folderMessage = 'Select directory with data to analyze'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase, folderMessage);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '.mat']);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
rateMaps = data; rateMapSettings = settings; 

% Settings: 
rateMapSettings.firingRates.meanThresh = 0.1; 
rateMapSettings.firingRates.maxThresh = 1; 
rateMapSettings.firingRates.plot.display = 'no'; 

% Outputs:
dataByFiringRate = splitLowAndHighFR_v01_20240802(rateMaps, rateMapSettings, processedDataFolder); 
plotFiringRateMetrics_v1_20250521(dataByFiringRate, rateMapSettings)

%% Step 9: Get spatial metrics and associated barcode (takes seconds)
%clear;clc;tic;

% Inputs:
fileNameBase = 'rateMapsByFiringRate';
folderMessage = 'Select directory with data to analyze'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase, folderMessage);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}(1:end-4)];
load([filePath{1}, '\', loadFileName, '.mat']);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
rateMapByFRStructure = data; rateMapByFRSettings = settings; 

% Settings: 
rateMapByFRSettings.rateMaps.lowThresh = 0.3; % fields will be counted as contiguous bins above 10% of max
rateMapByFRSettings.rateMaps.highThresh = 0.5; % fields must include one bin that is above 50% of the max
rateMapByFRSettings.rateMaps.biggerFieldModification = 2; % for use in phase precession control analysis: extend field by X bins
rateMapByFRSettings.rateMaps.smallerFieldModification = -2; % for use in phase precession control analysis: reduce field by X bins
rateMapByFRSettings.rateMaps.plot.display = 'yes'; 
rateMapByFRSettings.rateMaps.plot.genotypes = 'all'; % Could be 'WT', 'KO', or 'all'
rateMapByFRSettings.rateMaps.plot.animals = 'all'; % Could be a number between 1 and 8 or 'all'
rateMapByFRSettings.rateMaps.plot.cells = 'all'; % Could be a number between 1 and 8 or 'all'
rateMapByFRSettings.rateMaps.plot.direction = 'all'; % Could be 1, 2, or 'all'

% Outputs:
spatialMetrics = getSpatialMetrics_v1_20240724(rateMapByFRStructure, rateMapByFRSettings, processedDataFolder); toc
spatialMetrics = controlSpatialMetrics_v1_20250425(spatialMetrics, rateMapByFRSettings); 
plotSpatialMetrics_v1_20250305(spatialMetrics, rateMapByFRSettings); 
    
%% Step 10: Get the in-field spikes for each field (takes seconds)
clear;clc;tic;

% Inputs: 
fileNameBase = 'spatialMetrics';
folderMessage = 'Select directory with data to analyze'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase, folderMessage);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
binnedSpikes = data; binnedSpikeSettings = settings;

% Settings: 

% Outputs: 
inFieldSpkTimes = getInFieldSpikes_v1_20240805(binnedSpikes, binnedSpikeSettings, processedDataFolder); toc

%% Step 11: Additional rate map analysis

%% Step 12: Stability analysis

%% Step 13: Get the theta modulation (takes ~1.5 hours)
clear;clc;tic;
display(['Estimated time to finish is ', datestr(datetime('now')+hours(1.5), 'HH:MM:SS')]); 

% Inputs: 
fileNameBase = 'inFieldSpkTimes';
folderMessage = 'Select directory with data to analyze'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase, folderMessage);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
inFieldSpkTimes = data; thetaSettings = settings;

% Settings: 
thetaSettings.theta.frequencyBand = [4,12]; 
thetaSettings.theta.fieldsToAnalyze = 'best field'; 
thetaSettings.theta.numBins = 24; % Number of bins for rose plots
thetaSettings.theta.plot.display = 'no'; % Do you want to generate plots?
thetaSettings.theta.plot.genotypes = 2; % If plotting, for what genotype? Could be 'all'
thetaSettings.theta.plot.animals = 1; % If plotting, for what animals? Could be 'all'
thetaSettings.theta.plot.cells = 7; % If plotting, for what cells? Could be 'all'
thetaSettings.theta.plot.direction = 1; % If plotting, for what direction? Could be 'all'

% Outputs: 
thetaData = getThetaModulation_v1_20240806(inFieldSpkTimes, thetaSettings, processedDataFolder); toc
%controlsForTheta_v1_20250228(thetaData, thetaSettings); % Nothing in here
plotThetaModulation_v1_20241216(thetaData, thetaSettings); 

%% Step 14: Get the phase precession (takes ~8 min)
clear;clc;close all; tic;

% Inputs: 
fileNameBase = 'theta';
folderMessage = 'Select directory with data to analyze'; 
filePath = getMostRecentFilePath_v1_20240723(fileNameBase, folderMessage);
loadFileName = [fileNameBase, '_v', filePath{2}, filePath{3}];
load([filePath{1}, '\', loadFileName]);
[processedDataFolder, ~, ~] = fileparts(filePath{1});
thetaData = data; phasePrecessionSettings = settings;

% Settings: 
phasePrecessionSettings.phasePrecession.spatialBinThreshold = 0; % minimum number of spatial bins needed
phasePrecessionSettings.phasePrecession.spikeThreshold = 5; % minimum number of spikes needed
phasePrecessionSettings.phasePrecession.slopeRange = [-4:0.125:4]; % range of slopes to try to fit 
phasePrecessionSettings.phasePrecession.significanceThreshold = 0.05; % maximum acceptable significance level of line fit
phasePrecessionSettings.phasePrecession.trialThreshold = 3; % minimum number of trials
phasePrecessionSettings.phasePrecession.ISIthreshold = 1000; % max time between spikes in msec
phasePrecessionSettings.phasePrecession.fieldsToAnalyze = 'best field'; % Could also be 'best field'
phasePrecessionSettings.phasePrecession.positionType = 'unbinned';
phasePrecessionSettings.phasePrecession.normalized = 'no'; % Is the field normalized?
phasePrecessionSettings.phasePrecession.circularity = 'shift'; % How is circularity accounted for? Could also be set to 'shift'
phasePrecessionSettings.phasePrecession.fit = 'linearFit'; % Using a linear or circular fit?
phasePrecessionSettings.phasePrecession.timeRange = 3*125; % Over what range of time should spikes occur? 625 msec = 5 theta cycles
phasePrecessionSettings.phasePrecession.control.bootstrapN = 1000; % Number of iterations to bootstrap over
phasePrecessionSettings.phasePrecession.plot.display = 'yes'; % Determines which plots will be generated; 'yes' plots all while 'relationshipsOnly' only plots population data
phasePrecessionSettings.phasePrecession.plot.genotypes = 'all'; % If plotting, for what genotype? Could be 'all'
phasePrecessionSettings.phasePrecession.plot.animals = 'all'; % If plotting, for what animals? Could be 'all'
phasePrecessionSettings.phasePrecession.plot.cells = 'all'; % If plotting, for what cells? Could be 'all'
phasePrecessionSettings.phasePrecession.plot.direction = 'all'; % If plotting, for what direction? Could be 'all'

% Outputs: 
phasePrecessionData = getPhasePrecession_v1_20240806(thetaData, phasePrecessionSettings, processedDataFolder); toc 

%%
clc; close all; %phasePrecessionData = data;
phasePrecessionData = controlPhasePrecession_v1_20250228(phasePrecessionData, phasePrecessionSettings); 
%%
clc;
%phasePrecessionData = data; phasePrecessionSettings = settings; 
plotPhasePrecession_v2_20250305(phasePrecessionData, phasePrecessionSettings); 

%% Playing around
% Spearman's Correlation
%fieldSizes_wt = phasePrecessionData.populationData(1).phasePrecession.averageFieldSizes;
%fieldSizes_ko = phasePrecessionData.populationData(2).phasePrecession.averageFieldSizes;
fieldSizes_wt = phasePrecessionData.populationData(1).phasePrecession.MedianFieldSizes;
fieldSizes_ko = phasePrecessionData.populationData(2).phasePrecession.MedianFieldSizes;
slopes_wt = phasePrecessionData.populationData(1).phasePrecession.Slopes;
slopes_ko = phasePrecessionData.populationData(2).phasePrecession.Slopes;
validIdx = ~isnan(fieldSizes_wt) & ~isnan(slopes_wt);
field_sizes_wt = fieldSizes_wt(validIdx);
phase_slopes_wt = slopes_wt(validIdx);
validIdx = ~isnan(fieldSizes_ko) & ~isnan(slopes_ko);
field_sizes_ko = fieldSizes_ko(validIdx);
phase_slopes_ko = slopes_ko(validIdx);
format long
[rho, pval] = corr([(field_sizes_wt)';(field_sizes_ko)'], [phase_slopes_wt'; phase_slopes_ko'], 'Type', 'Spearman')
fprintf('Spearman correlation p = %.10f\n', pval);
[rho, pval] = corrcoef([log(field_sizes_wt)';log(field_sizes_ko)'], [phase_slopes_wt'; phase_slopes_ko']);
fprintf('Pearsons correlation p = %.10f\n', pval(2));

% Results suggest strongly correlated as we might expect

% Does log transform create linear relationship? 
figure;
scatter(log(field_sizes_wt), phase_slopes_wt)
scatter(log(field_sizes_ko), phase_slopes_ko, 'g')
% It looks fairly linear

% Log-transformed regression
wt_genotype = ones(1,length(field_sizes_wt)); 
ko_genotype = 2*ones(1,length(field_sizes_ko)); 
Genotype = [wt_genotype, ko_genotype]; 
log_FieldSize = [log(field_sizes_wt), log(field_sizes_ko)];
phase_slopes = [phase_slopes_wt, phase_slopes_ko]; 
tbl = table(phase_slopes', log_FieldSize', Genotype', 'VariableNames', {'Slope', 'LogFieldSize', 'Genotype'});
mdl = fitlm(tbl, 'Slope ~ LogFieldSize + Genotype');
disp(mdl)
% The logFieldSize p-value is highly significant suggesting that field size
% plays an important role in the slope
% The genotype is not significant suggesting that genotype does not play an
% important role in the slope

%%
% Same log-transformed regression but looking for an interaction between
% gentoype and size
mdl = fitlm(tbl, 'Slope ~ LogFieldSize * Genotype');
disp(mdl)
% This model only explains 22% of the variance which is lower than the
% prevous model so we should stick with the previous model

% Now trying to tie in theta modulation
% Is theta modulation related to phase precession? 
theta_wt = phasePrecessionData.populationData(1).MVL; 
theta_ko = phasePrecessionData.populationData(2).MVL; 
figure; scatter(theta_wt, slopes_wt); hold on; scatter(theta_ko, slopes_ko); 
xlabel('theta'); ylabel('slopes');
% It doesn't look like there's a relationship

% Is field size related to theta
figure; scatter(theta_wt, fieldSizes_wt); hold on; scatter(theta_ko, fieldSizes_ko); 

% Include it in the regression
validIdx = ~isnan(fieldSizes_wt) & ~isnan(slopes_wt);
theta_wt_noNaN = theta_wt(validIdx);
validIdx = ~isnan(fieldSizes_ko) & ~isnan(slopes_ko);
theta_ko_noNaN = theta_ko(validIdx);
theta = [theta_wt_noNaN, theta_ko_noNaN]; 

[rho, pval] = corr(theta', [(field_sizes_wt)';(field_sizes_ko)'], 'Type', 'Spearman')

%tbl = table(phase_slopes', log_FieldSize', theta', Genotype', ...
%    'VariableNames', {'Slope', 'LogFieldSize', 'ThetaModulation', 'Genotype'});
tbl = table(phase_slopes', theta', Genotype', ...
    'VariableNames', {'Slope', 'ThetaModulation', 'Genotype'});
mdl = fitlm(tbl, 'Slope ~ ThetaModulation + Genotype');
disp(mdl)
%{
tbl = table(phase_slopes', theta', Genotype', ...
    'VariableNames', {'Slope', 'ThetaModulation', 'Genotype'});
mdl = fitlm(tbl, 'Slope ~ ThetaModulation + Genotype');
disp(mdl)
%}

%%
% Now I'm curious about bursting
% Is bursting related to theta modulation? 
dataForBursts = phasePrecessionData.cellData;
genotypes = fieldnames(dataForBursts);
for iGenotype = 1:length(genotypes); 
    genotypeData = dataForBursts.(genotypes{iGenotype}).highFiring; 
    burstIndex = []; 
    for iAnimal = 1:length(genotypeData)
        for iCluster = 1:length(genotypeData{iAnimal})
            if ~isempty(genotypeData{iAnimal}(iCluster).spikeTiming)
                bursts = genotypeData{iAnimal}(iCluster).spikeTiming.bursts;
                directions = fieldnames(bursts); 
                for iDir = 1:length(directions); 
                    burstsByDir = bursts.(directions{iDir}); 
                    allBursts = []; 
                    for iTrial = 1:length(burstsByDir); 
                        allBursts = [allBursts, burstsByDir{iTrial}]; 
                    end
                    burstIndex = [burstIndex, sum(~isnan(allBursts))/length(allBursts)]; 
                end
            end
        end
    end
    if iGenotype == 1; 
        burstIndex_wt = burstIndex; 
    elseif iGenotype == 2;
        burstIndex_ko = burstIndex; 
    end
end
figure; scatter(burstIndex_wt, theta_wt); hold on; scatter(burstIndex_ko, theta_ko);
xlabel('Bursting'); ylabel('Theta');  

% Is bursting related to phase precession? 
figure; scatter(burstIndex_wt, slopes_wt); hold on; scatter(burstIndex_ko, slopes_ko);
xlabel('Bursting'); ylabel('Slopes');  

% Is theta moudlation related to phase precession? 
figure; scatter(theta_wt, slopes_wt); hold on; scatter(theta_ko, slopes_ko);
xlabel('Theta'); ylabel('Slopes');  

% A big ole linear regression model
validIdx = ~isnan(fieldSizes_wt) & ~isnan(slopes_wt);
bursting_wt_noNaN = burstIndex_wt(validIdx);
validIdx = ~isnan(fieldSizes_ko) & ~isnan(slopes_ko);
bursting_ko_noNaN = burstIndex_ko(validIdx);
bursting = [bursting_wt_noNaN, bursting_ko_noNaN];

[rho, pval] = corr(bursting', phase_slopes', 'Type', 'Spearman');
[rho, pval] = corr(theta', phase_slopes', 'Type', 'Spearman')
[rho, pval] = corr(theta', bursting', 'Type', 'Spearman');
[h, p] = corrcoef(theta, phase_slopes)


%{
tbl = table(phase_slopes', log_FieldSize', bursting', Genotype', ...
    'VariableNames', {'Slope', 'FieldSize', 'Bursting', 'Genotype'});
mdl = fitlm(tbl, 'Slope ~ FieldSize + Bursting + Genotype');
disp(mdl)
%}
%{
tbl = table(phase_slopes', theta', log_FieldSize', bursting', Genotype', ...
    'VariableNames', {'Slope', 'Theta', 'FieldSize', 'Bursting', 'Genotype'});
mdl = fitlm(tbl, 'Slope ~ Theta + FieldSize + Bursting + Genotype');
disp(mdl)
%}
tbl = table(phase_slopes', theta', log_FieldSize', Genotype', ...
    'VariableNames', {'Slope', 'Theta', 'FieldSize', 'Genotype'});
mdl = fitlm(tbl, 'Slope ~ Theta + FieldSize + Genotype');
disp(mdl)
R2 = mdl.Rsquared.Ordinary;
vif_fieldSize = 1 / (1 - R2)

%%
% Is bursting related to theta modulation in a model? 
tbl = table(theta', bursting', Genotype', ...
    'VariableNames', {'ThetaModulation', 'Bursting', 'Genotype'});
mdl = fitlm(tbl, 'ThetaModulation ~ Bursting + Genotype');
disp(mdl)
%% 
% Does genotype alone contribute to phase precession slope? 
tbl = table(phase_slopes', Genotype', ...
    'VariableNames', {'Slope', 'Genotype'});
mdl = fitlm(tbl, 'Slope ~ Genotype');
disp(mdl)

%% Using a model with all the interactions
tbl = table(phase_slopes', Genotype', log_FieldSize', theta', ...
    'VariableNames', {'Slope', 'Genotype', 'FieldSize', 'ThetaMod'});

mdl = fitlm(tbl, 'Slope ~ Genotype * FieldSize + Genotype * ThetaMod');
disp(mdl)

%% genotype only
slopes_wt = phasePrecessionData.populationData(1).phasePrecession.Slopes;
slopes_ko = phasePrecessionData.populationData(2).phasePrecession.Slopes;
phase_slopes = [slopes_wt, slopes_ko]; 
wt_genotype = ones(1,length(slopes_wt)); 
ko_genotype = 2*ones(1,length(slopes_ko)); 
Genotype = [wt_genotype, ko_genotype]; 
validIdx = ~isnan(phase_slopes) & ~isnan(Genotype);
phase_slopes = phase_slopes(validIdx); 
Genotype = Genotype(validIdx); 

tbl = table(phase_slopes', Genotype', ...
    'VariableNames', {'Slope', 'Genotype'});

mdl = fitlm(tbl, 'Slope ~ Genotype');
disp(mdl)

%% field size only
fieldSizes_wt = phasePrecessionData.populationData(1).phasePrecession.MedianFieldSizes;
fieldSizes_ko = phasePrecessionData.populationData(2).phasePrecession.MedianFieldSizes;
slopes_wt = phasePrecessionData.populationData(1).phasePrecession.Slopes;
slopes_ko = phasePrecessionData.populationData(2).phasePrecession.Slopes;
validIdx = ~isnan(fieldSizes_wt) & ~isnan(slopes_wt);
field_sizes_wt = fieldSizes_wt(validIdx);
phase_slopes_wt = slopes_wt(validIdx);
validIdx = ~isnan(fieldSizes_ko) & ~isnan(slopes_ko);
field_sizes_ko = fieldSizes_ko(validIdx);
phase_slopes_ko = slopes_ko(validIdx);

% Log-transformed regression
log_FieldSize = [log(field_sizes_wt), log(field_sizes_ko)];
phase_slopes = [phase_slopes_wt, phase_slopes_ko]; 
tbl = table(phase_slopes', log_FieldSize', 'VariableNames', {'Slope', 'LogFieldSize'});
mdl = fitlm(tbl, 'Slope ~ LogFieldSize');
disp(mdl)

%% theta mvl only
theta_wt = phasePrecessionData.populationData(1).MVL; 
theta_ko = phasePrecessionData.populationData(2).MVL;
slopes_wt = phasePrecessionData.populationData(1).phasePrecession.Slopes;
slopes_ko = phasePrecessionData.populationData(2).phasePrecession.Slopes;
validIdx = ~isnan(theta_wt) & ~isnan(slopes_wt);
theta_wt = theta_wt(validIdx);
phase_slopes_wt = slopes_wt(validIdx);
validIdx = ~isnan(theta_ko) & ~isnan(slopes_ko);
theta_ko = theta_ko(validIdx);
phase_slopes_ko = slopes_ko(validIdx);

mvl = [theta_wt, theta_ko]; 
phase_slopes = [phase_slopes_wt, phase_slopes_ko]; 
tbl = table(phase_slopes', mvl', 'VariableNames', {'Slope', 'MVL'});
mdl = fitlm(tbl, 'Slope ~ MVL');
disp(mdl)

%% genotype + field size + mvl
slopes_wt = phasePrecessionData.populationData(1).phasePrecession.Slopes;
slopes_ko = phasePrecessionData.populationData(2).phasePrecession.Slopes;
phase_slopes = [slopes_wt, slopes_ko]; 
wt_genotype = ones(1,length(slopes_wt)); 
ko_genotype = 2*ones(1,length(slopes_ko)); 
Genotype = [wt_genotype, ko_genotype]; 
fieldSizes_wt = phasePrecessionData.populationData(1).phasePrecession.MedianFieldSizes;
fieldSizes_ko = phasePrecessionData.populationData(2).phasePrecession.MedianFieldSizes;
fieldSizes = [fieldSizes_wt,fieldSizes_ko]; 
theta_wt = phasePrecessionData.populationData(1).MVL; 
theta_ko = phasePrecessionData.populationData(2).MVL;
theta = [theta_wt, theta_ko]; 

validIdx = ~isnan(phase_slopes) & ~isnan(Genotype) & ~isnan(fieldSizes) & ~isnan(theta);
phase_slopes = phase_slopes(validIdx); 
Genotype = Genotype(validIdx); 
log_FieldSize = log(fieldSizes(validIdx)); 
theta = theta(validIdx); 

tbl = table(phase_slopes', Genotype', log_FieldSize', theta', ...
    'VariableNames', {'Slope', 'Genotype', 'FieldSize', 'Theta'});
mdl = fitlm(tbl, 'Slope ~ Genotype + FieldSize + Theta');
disp(mdl)

% VIFs for full model
mdl_fs = fitlm(tbl, 'FieldSize ~ Genotype + Theta');
R2_fs = mdl_fs.Rsquared.Ordinary;
vif_FieldSize = 1 / (1 - R2_fs)

mdl_theta = fitlm(tbl, 'Theta ~ Genotype + FieldSize');
R2_theta = mdl_theta.Rsquared.Ordinary;
vif_Theta = 1 / (1 - R2_theta)

mdl_geno = fitlm(tbl, 'Genotype ~ FieldSize + Theta');
R2_geno = mdl_geno.Rsquared.Ordinary;
vif_Genotype = 1 / (1 - R2_geno)

%% Interaction Full Model
tbl = table(phase_slopes', Genotype', log_FieldSize', theta', ...
    'VariableNames', {'Slope', 'Genotype', 'FieldSize', 'ThetaMod'});

mdl = fitlm(tbl, 'Slope ~ Genotype * FieldSize + Genotype * ThetaMod');
disp(mdl)

