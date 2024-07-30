function data = calculateRateMap_v1_20240718(data, settings)
    % Calculate the linearized rate maps across all trials for each cells
    % Written by Anja Payne
    % Last Modified: 07/29/2024

    % Inputs:
    %   1) data: the matlab structure where the spikes (binned into 4cm
    %      bins and separated into trials/direction) are stored
    %   2) settings: the matlab settings file where the ______ thresholds 
    %      are stored
    
    % Outputs:
    %   1) data: the structure with the high-velocity spike times included

    % Steps:
    %   1) Loop through genotypes, animals, and cells and calculate the 
    %      rate maps for each trial
    %   2) If the user has selected this option, plot the rate maps as a
    %      two-dimensional heat map with time on the x-axis and trial
    %      number on the y-axis. 
    %   3) Ask the user if they want to save the newly generated data
    
    nBins = settings.rateMaps.trackSize / settings.rateMaps.binSize;
    for iGenotype = 1:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype}); 
        for iAnimal = 1:length(genotypeData); 
            if isempty(genotypeData{iAnimal}) == 1; 
                continue
            else
                [~,n] = size(genotypeData{iAnimal});
                for iCluster = 1:n;
                    %% Step 1: Get the rate maps
                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                    for iDir = 1:2;
                        % Assign binned positions and spikes
                        if iDir == 1; 
                            posBinned = genotypeData{iAnimal}(iCluster).binnedPosByTrial.cw;
                            spikesBinned = genotypeData{iAnimal}(iCluster).binnedSpikesByTrial.cw;
                        elseif iDir == 2; 
                            posBinned = genotypeData{iAnimal}(iCluster).binnedPosByTrial.ccw;
                            spikesBinned = genotypeData{iAnimal}(iCluster).binnedSpikesByTrial.ccw;
                        end
                        map = []; timeMap = []; 
                        for iTrials = 1:length(spikesBinned); 
                            [map(iTrials,:), timeMap(iTrials,:), spikeMap(iTrials,:)] = calculateLinearizedRateMap(posBinned{iTrials}, spikesBinned{iTrials}, nBins);
                        end
                        % Shift the map so that the reward zone is at the end of the track
                        map = circshift(map, -2, 2);
                        % Multiple the map by the sampling rate to get
                        % spikes per second
                        map = map * settings.velocity.samplingRate; 
                        
                        % Append to data structure
                        if iDir == 1; 
                            map = fliplr(map); % so that running proceeds left to right
                            timeMap = fliplr(timeMap); 
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.rateMap.cw = map;
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.timeMap.cw = timeMap;
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.trialAverageRates.cw = nanmean(map,1);
                        elseif iDir == 2; 
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.rateMap.ccw = map;
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.timeMap.ccw = timeMap;
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.trialAverageRates.ccw = nanmean(map,1);
                        end
                    end
                end
            end
        end
    end
    
    %% Step 2: Step 2: Save
    saveFile_v1_20240718(data, settings, 'rateMaps') 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [map, timeMap, spikeMap] = calculateLinearizedRateMap(pos, spikes, nBins)
    % Gets the rate map for a linearized representation of the map
    % Original code from Leutgeb lab with modifications by A. Payne
    % Inputs:
    %   1) pos: the binned positions
    %   2) spikes: the binned spikes
    % Outputs:
    %   1) map: a heatmap of firing rate per bin
    %   2) timeMap: a heatmap of time per bin
    %   3) spikeMap: a heatmap of spike number per bin
    
    spikeMap = zeros(1, nBins); 
    for iSpike = 1:length(spikes);
        spikeMap(spikes(iSpike)) = spikeMap(spikes(iSpike)) + 1; 
    end
    timeMap = zeros(1, nBins); 
    for iTime = 1:length(pos);
       timeMap(pos(iTime)) = timeMap(pos(iTime)) + 1; 
    end

    % Set unoccupied pixels to NaN; 
    minOccupancy = 1; 
    timeMap(timeMap < minOccupancy) = NaN; 
    rawMap = spikeMap./timeMap; 
    % Shift the map so that the reward zone is at the end
    rawMap = circshift(rawMap, -2, 2);
    % Add some padding for smoothing to work effectively
    rawMap = padarray(rawMap', 2, NaN); 
    
    % Smooth the rate map using the Gaussian filter as well as an adaptive
    % smoothing technique
    % Define the Gaussian smoothing that will be applied to each bin. 
    box = [0.0200 0.1000 0.1600 0.1000 0.0200]; 
    for i = 3:length(rawMap)-2; 
        current_sum = 0;
        current_box = 0;
        box_i = 0;
        for ii = i-2:i+2;
          box_i = box_i + 1;
          if ~isnan(rawMap(ii));
             current_sum = current_sum + rawMap(ii) * box(box_i);
             current_box = current_box + box(box_i); 
          end
        end
        map(i) = current_sum/current_box; 
    end
    
    map = map(3:end); 
    
end

