function data = getStability_v1_20250608(data, settings, processedDataPath)
    % Gets the stability and stability-related analysis from the rate maps
    % Written by Anja Payne
    % Last Modified: 06/08/2026
    % Inputs:
    %   1) data: the matlab structure where the trial-averaged rate
    %      maps are stored
    %   2) settings: the matlab settings file where the thresholds are
    %      stored
    
    % Outputs:
    %   1) data: the structure with stability metrics appended

    % Steps:
    %   1) Loop through genotypes, animals, and cells and get the
    %   correlation between epochs for the full track
    %   2) Get the correlation between epochs for the in-field bins
    %   3) Get the correlation between epochs for the out-of-field bins
    
    %%%%%%%%
    
    % Note for myself: For stability, I ran
    % Z:\Anja\Matlab Code\CreatePlots\Workflows\Stability\AP_Workflow_PlotSequentialEpochStability_BothDir_Sparse.m
    % for the full track. 
    % For the full track shuffle I modified the function referenced in that
    % workflow (splitEpochs_GetSequentialCorr_CollapsedMaps)
    
    
    
    
    %% Step 1: Get the correlation between epochs for the full track

    for iGenotype = 1:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype});
        for iFR = 1:length(fieldnames(genotypeData)); 
            FRoptions = fieldnames(genotypeData); 
            FRdata = genotypeData.(FRoptions{iFR}); 
            for iAnimal = 1:length(FRdata); 
                if isempty(FRdata{iAnimal}) == 1; 
                    continue
                else
                    [~,n] = size(FRdata{iAnimal});
                    for iCluster = 1:n;
                        if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                            continue
                        else
                            display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                            directions = fieldnames(FRdata{iAnimal}(iCluster).rateMaps.rateMap);
                            for iDir = 1:length(directions);
                                % Load the file with the binned position
                                load(genotypeData{iAnimal}(iCluster).posBinFile); 
                                
                                
                                % Extract data based on running direction
                                outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir), 'stability');
                                map = outputData.map; timeMap = outputData.timeMap; 
                                trialAveragedMap = nanmean(map, 1); 
                                trialAveragedTimeMap = nanmean(timeMap, 1); 
                            end
                        end
                    end
                end
            end
        end
    end
end

                                
                                
                                
                                
                                
                                
                                