function calculateRateMap_v1_20240718
% Calculate the linearized rate maps across all trials for each cells
% Written by Anja Payne
% Last Modified: 07/18/2024

% Inputs:

% Outputs:

% Steps:
%   1) Load relevant data and define necessary pathways
%   2) Loop through each animal and each cell. Center the x and y position 
%      around 0 and convert coordinates to mm
%   3) Exclude spikes that occur before the video, after the video, and
%      during periods of low velocity. 
%   4) Split the file into trials
%   5) Linearize the track into bins and assign each spike to the nearest
%      bin. 
%   6) Calculate the rate maps for each trial
%   7) Plot as a 2 dimensional heat map with time on the x-axis and trial
%      number on the y-axis. 


%% Step 1: %%%%
clear
clc
% Define pathways
%addpath(genpath('Z:\Anja\Matlab Code\CreateRateMaps')); 
%addpath(genpath('Z:\Anja\Matlab Code\LoadFiles')); 
%addpath(genpath('Z:\Anja\Matlab Code\GeoffCode\Organized Code'));
animalFile = 'AP_AnimalData_Cohort13_Track'; 
animalInfoFile = ['Z:\Anja\Data\In Vivo Data\AnimalInfo_ExcelFiles\', animalFile]; 
animalInfo = excel2Mat(animalInfoFile)
animalInfo_sessionSpecific = selectSessionType(animalInfo, 'Track'); 
%%
% Define parameters
trackWidth = 52; % In cm
trackLength = 80; % In cm
close all; 
%%%% Step 2: %%%%
for iAnimal = 1%:length(animalInfo_sessionSpecific); 
   for iCluster = 64%1:length(animalInfo_sessionSpecific(iAnimal).directory);
       display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
       %%%% Step 1: %%%%
        % Load the event file
        event_file = [animalInfo_sessionSpecific(iAnimal).directory{iCluster}, '\Events.nev'];
        % Load the spike data for the relevant cluster
        load([animalInfo_sessionSpecific(iAnimal).directory{iCluster}, '\', animalInfo_sessionSpecific(iAnimal).clusters{iCluster}]);
        % Load the video file
        idcs = strfind(animalInfo_sessionSpecific(iAnimal).directory{iCluster}, '\');
        parent_directory = animalInfo_sessionSpecific(iAnimal).directory{iCluster}(1:idcs(end)-1);
        load([parent_directory, '\Position.mat']);
        % Load the csc data
        ephys_directory = animalInfo_sessionSpecific(iAnimal).directory{iCluster}; 
        CSC_file = [ephys_directory, '\CSC1.ncs']; 
        csc_data = readCSC(CSC_file); 
        csc_timepoints_lowSampling = csc_data.ts; 
        csc_samples = csc_data.samp; 
        
        %%%% Step 2: %%%%
        % Exclude any position data that occurs before the video turns on. 
        clipped_indices = find(aligned_pos(:,1) ~= -1000);
        x_clipped = aligned_pos(clipped_indices, 2);
        y_clipped = aligned_pos(clipped_indices, 1); 
        t_clipped = aligned_pos_ts(clipped_indices); 
        % Center the position data around the coordinates (0,0)
        [x_pos_centered, y_pos_centered] = centerBox(x_clipped, y_clipped); 
        % Convert the scale from pixels to cm
        widthConversionFactor = trackWidth/(max(x_pos_centered) - min(x_pos_centered)); 
        lengthConversionFactor = trackLength/(max(y_pos_centered) - min(y_pos_centered)); 
        x_pos_cm = x_pos_centered*widthConversionFactor;
        y_pos_cm = y_pos_centered*lengthConversionFactor; 
        
        %%%% Step 3: %%%%
        % Interpolate the CSC timepoints so that they are at the same
        % sampling rate as the CSC samples
        csc_samples_points = linspace(1, length(csc_timepoints_lowSampling), length(csc_samples)); 
        csc_timepoints_lowSampling_points = [1:length(csc_timepoints_lowSampling)]; 
        csc_timepoints = interp1(csc_timepoints_lowSampling_points, csc_timepoints_lowSampling, csc_samples_points, 'linear'); 
        csc_timepoints_msec = csc_timepoints/1000; 
        ephys_t0 = csc_timepoints_msec(1); 
        % Exclude any spikes that occur before or after the video has started. 
        spike_time_msec = TS/10; 
        %exclude_timepoints = csc_timepoints_msec(binaryVideoOn == 0); 
        % Exclude any spikes that occur during periods of low velocity. 
        [include_spikes, ave_vel, std_vel, vel_thresh, velocity] = excludeLowVelocity(x_pos_cm, y_pos_cm, t_clipped, spike_time_msec, ephys_t0); 

        %%%% Step 4: %%%%
        % Split the spikes and position into trials. 
        [xTrials,yTrials,tTrials,sTrials] = splitFilesByTrials(include_spikes, x_pos_cm, y_pos_cm, t_clipped);

        %%%% Step 5: %%%%
        % Determine how you want the linearized trials organized. A value
        % of 0 will organize the trials in time (1 to the end). A value of
        % 1 will organize the trials into clockwise vs. counterclockwise. 
        organization = 1;
        % Linearize the track into bins and assign each spike to the
        % nearest bin
        trials = {};
        [occupancyBinned, spikesBinned, trials] = linearizeData(xTrials, yTrials, tTrials, sTrials, organization); 
%%
        %%%% Step 6: %%%%
        % Calculate the rate maps
        map = []; timeMap = []; 
        for iTrials = 1:length(spikesBinned); 
            [map(iTrials,:), timeMap(iTrials,:), spikeMap(iTrials,:)] = calcRateMapLinearized(occupancyBinned{iTrials}, spikesBinned{iTrials});
        end
        % Shift the map so that the reward zone is at the end of the track
        %map = circshift(temp_map, -2, 2);
        
        %%
        %%%% Step 7: %%%%
        % Plot as a 2 dimensional heat map with time on the x-axis and trial
        % number on the y-axis.
        %map = map(1:length(trials)/2,:);
        imAlpha = ones(size(map)); 
        imAlpha(isnan(map)) = 0;
        h1 = figure(1);
        imagesc(map*16000, 'AlphaData', imAlpha);
        set(gca, 'color', [0.8,0.8,0.8]); 
        title(['Linearized Rate Map for ', animalInfo_sessionSpecific(iAnimal).cell_name{iCluster}], 'Interpreter', 'none');
        set(gca, 'ytick', [1:length(trials)]); 
        set(gca, 'YTickLabels', trials); 
        xlabel('Bin');
        ylabel('Trial'); 
        colorbar;
        %saveas(h1, ['Z:\Anja\Figures\4. Spatial Characteristics\LinearTrack\Linearized\Sparse\PotentialExamples\', animalInfo_sessionSpecific(iAnimal).cell_name{iCluster}, '_AverageLinearizedMap.tif'],'tif'); 
        %saveas(h1, [parent_directory, '\', animalInfo_sessionSpecific(iAnimal).cell_name{iCluster}, '_BinSize4_LinearizedBurstRateMap.fig'],'fig');
        %saveas(h1, ['Z:\Anja\Data\Cells\Rate Maps\Linearized Rate Maps\Linear Track\', animalInfo_sessionSpecific(iAnimal).cell_name{iCluster}, '_BinSize4_Linearized.tif'],'tif');       
        %clf
        % If desired, plot the average rate map
        
        h2 = figure(2);
        h2 = imagesc(nanmean(map(1:length(trials)/2,:))*16000);
        %saveas(h2, ['Z:\Anja\Figures\4. Spatial Characteristics\LinearTrack\Linearized\Dense\PotentialExamples\bestKOexamples\', animalInfo_sessionSpecific(iAnimal).cell_name{iCluster}, '_AverageLinearizedMap.tif'],'tif'); 
        %%
        
        % Concatenate all rate maps into one variable
        allRateMaps{iAnimal, iCluster} = map; 
        % Concatenate all trials into one variable
        allTrials{iAnimal, iCluster} = trials; 
        % Concatenate all time maps into one variable
        allTimeMaps{iAnimal, iCluster} = timeMap; 

   end
end

% Save variables
% idcs   = strfind(animalInfoFile,'\');
% save_directory = animalInfoFile(idcs(end)+1:end);
% save(['Z:\Anja\Data\In Vivo Data\AnalyzedData\', save_directory, 'Velocity2_BinSize4_LinearizedRateMaps'], 'allRateMaps', 'allTrials', 'allTimeMaps');
%save(['Z:\Anja\Data\In Vivo Data\AnalyzedData\', save_directory, '_V2B4_Bursts'], 'allRateMaps', 'allTrials', 'allTimeMaps');
%save(['Z:\Anja\Data\In Vivo Data\AnalyzedData\', save_directory, '_V2B4_Events'], 'allRateMaps', 'allTrials', 'allTimeMaps');
%save(['Z:\Anja\Data\In Vivo Data\AnalyzedData\', save_directory, '_V2B4_Singles'], 'allRateMaps', 'allTrials', 'allTimeMaps');

%% Temporary plotting
trial = 11; trial2 = 41; 
close all; figure(1); subplot(1,2,1); hold on;
plot(xTrials{trial}, yTrials{trial}, 'k'); 
plot(xTrials{trial}(1:100000), yTrials{trial}(1:100000), 'g'); 
[~, ~, sPoints] = intersect(sTrials{trial}, round(tTrials{trial}))
scatter(xTrials{trial}(sPoints), yTrials{trial}(sPoints), '*r'); 
subplot(1,2,2); 
imagesc(map(trial2,:)); 
%set(gca, 'ytick', [1:length(trials)]); 
set(gca, 'YTickLabels', trials{trial2}); 
set(gcf, 'position', [200, 200, 1000, 500]); 
%%
figure(1); clf; hold on;
plot([0:4:264], nanmean(map(length(trials)/2:end,:))*16000);
plot([1,264], [0.1*max(nanmean(map(length(trials)/2:end,:))*16000), 0.1*max(nanmean(map(length(trials)/2:end,:))*16000)]); 
plot([1,264], [0.5*max(nanmean(map(length(trials)/2:end,:))*16000), 0.5*max(nanmean(map(length(trials)/2:end,:))*16000)]); 