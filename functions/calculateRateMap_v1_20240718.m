function data = calculateRateMap_v1_20240718(data, settings, processedDataPath)
    % Calculate the linearized rate maps across all trials for each cells
    % Written by Anja Payne
    % Last Modified: 07/29/2024

    % Inputs:
    %   1) data: the matlab structure where the spikes (binned into 4cm
    %      bins and separated into trials/direction) are stored
    %   2) settings: the matlab settings file where the rate map thresholds 
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
                    % Load the position data
                    load(genotypeData{iAnimal}(iCluster).posBinFile); 
                    for iDir = 1:2;
                        % Assign binned positions and spikes
                        if iDir == 1; 
                            posBinned = binnedPosition.cw; 
                            spikesBinned = genotypeData{iAnimal}(iCluster).highVelocityData.spikePosBins.cw;
                        elseif iDir == 2; 
                            posBinned = binnedPosition.ccw;
                            spikesBinned = genotypeData{iAnimal}(iCluster).highVelocityData.spikePosBins.ccw;
                        end
                        map = []; timeMap = []; 
                        for iTrials = 1:length(spikesBinned); 
                            [map(iTrials,:), timeMap(iTrials,:), spikeMap(iTrials,:)] = calculateLinearizedRateMap(posBinned{iTrials}, spikesBinned{iTrials}, nBins);
                        end
                        clear posBinned;

                        % Multiple the map by the sampling rate to get
                        % spikes per second
                        map = map * settings.velocity.samplingRate; 
                        map = circshift(map, 1, 2); 
                        timeMap = timeMap / settings.velocity.samplingRate; 
                        timeMap = circshift(timeMap, 1, 2); 
                        
                        % Append to data structure
                        if iDir == 1; 
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.rateMap.cw = map;
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.timeMap.cw = timeMap;
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.trialAverageRates.cw = nanmean(map,1);
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).firingRates.max.cw = nanmax(nanmean(map,1)); 
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).firingRates.mean.cw = nanmean(nanmean(map,1)); 
                        elseif iDir == 2; 
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.rateMap.ccw = map;
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.timeMap.ccw = timeMap;
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).rateMap.trialAverageRates.ccw = nanmean(map,1);
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).firingRates.max.ccw = nanmax(nanmean(map,1)); 
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).firingRates.mean.ccw = nanmean(nanmean(map,1)); 
                        end
                    end
                end
            end
        end
    end
    
    %% Step 2: Step 2: Save
    saveFile_v1_20240718(processedDataPath, data, settings, 'rateMaps') 
    
end
