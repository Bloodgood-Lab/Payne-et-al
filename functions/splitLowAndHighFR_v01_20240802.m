function outputData = splitLowAndHighFR_v01_20240802(data, settings)
    % Given the spatial firing rates, split the data into low and high
    % firing rate cells
    % Written by Anja Payne
    % Last Modified: 08/02/2024

    % Inputs:
    %   1) data: the matlab structure where the rate maps, binned data, and
    %      firing rates are stored
    %   2) settings: the matlab settings file where the firing rate
    %      thresholds are stored
    
    % Outputs:
    %   1) data: the structure reorganized so that it is also split by low
    %      and high firing rates
    
    % Steps:
    %   1) Based on the firing rates, sort into low and high firing cells
    %   2) Save the data
    
    %% Step 1: Sort into low and high firing rate cells
    meanThresh = settings.firingRates.meanThresh;
    maxThresh = settings.firingRates.maxThresh; 
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
                            dir = 'cw';
                            meanFR = genotypeData{iAnimal}(iCluster).firingRates.mean.cw;
                            maxFR = genotypeData{iAnimal}(iCluster).firingRates.max.cw;
                            if meanFR > meanThresh || maxFR > maxThresh; 
                                fr = 'highFiring'; 
                            else 
                                fr = 'lowFiring'; 
                            end
                        elseif iDir == 2; 
                            dir = 'ccw';
                            meanFR = genotypeData{iAnimal}(iCluster).firingRates.mean.ccw;
                            maxFR = genotypeData{iAnimal}(iCluster).firingRates.max.ccw;
                            if meanFR > meanThresh || maxFR > maxThresh; 
                                fr = 'highFiring'; 
                            else 
                                fr = 'lowFiring'; 
                            end
                        end
                        
                        % Assign variables
                        outputData.(genotypes{iGenotype}).(fr){iAnimal}(iCluster).metaData.directory.(dir) = genotypeData{iAnimal}(iCluster).directory;
                        outputData.(genotypes{iGenotype}).(fr){iAnimal}(iCluster).metaData.fileName.(dir) = genotypeData{iAnimal}(iCluster).fileName;
                        outputData.(genotypes{iGenotype}).(fr){iAnimal}(iCluster).binnedSpikesByTrial.trials.(dir) = genotypeData{iAnimal}(iCluster).trials.(dir);
                        outputData.(genotypes{iGenotype}).(fr){iAnimal}(iCluster).binnedSpikesByTrial.binnedSpk.(dir) = genotypeData{iAnimal}(iCluster).spikePosBins.(dir);
                        outputData.(genotypes{iGenotype}).(fr){iAnimal}(iCluster).binnedSpikesByTrial.binnedPosFile = genotypeData{iAnimal}(iCluster).posBinFile;
                        outputData.(genotypes{iGenotype}).(fr){iAnimal}(iCluster).rateMaps.rateMap.(dir) = genotypeData{iAnimal}(iCluster).rateMap.rateMap.(dir);
                        outputData.(genotypes{iGenotype}).(fr){iAnimal}(iCluster).rateMaps.timeMap.(dir) = genotypeData{iAnimal}(iCluster).rateMap.timeMap.(dir);
                        outputData.(genotypes{iGenotype}).(fr){iAnimal}(iCluster).rateMaps.trialAverageMap.(dir) = genotypeData{iAnimal}(iCluster).rateMap.trialAverageRates.(dir);
                        outputData.(genotypes{iGenotype}).(fr){iAnimal}(iCluster).firingRates.mean.(dir) = genotypeData{iAnimal}(iCluster).firingRates.mean.(dir);
                        outputData.(genotypes{iGenotype}).(fr){iAnimal}(iCluster).firingRates.max.(dir) = genotypeData{iAnimal}(iCluster).firingRates.max.(dir);
                    end
                end
            end
        end
    end
    
    %% Step 2: Step 2: Save
    saveFile_v1_20240718(outputData, settings, 'rateMapsByFiringRate') 
    
end
