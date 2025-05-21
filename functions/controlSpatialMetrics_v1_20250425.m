function outputData = controlSpatialMetrics_v1_20250425(inputData, settings)
    % Performs control analysis related to the spatial metrics
    % Written by Anja Payne
    % Last Modified: 04/25/2025
    
    % Inputs:
    %   1) inputData: the matlab structure where the spatial metrics and
    %      rate maps have been stored
    %   2) settings: the settings file where the settings for the phase
    %      precession calculations are saved

    % Outputs:
    %   1) outputData: the modified matlab structure that was previously
    %      inputData but now with any control data included
    %     
    
    % Steps/Figures:
    %   1) 
    
    
    %% Step 1: 
    thresholds(1) = settings.rateMaps.lowThresh; 
    thresholds(2) = settings.rateMaps.highThresh;
    outputData = inputData; 
    genotypes = fieldnames(inputData);
    for iGenotype = 1:length(fieldnames(inputData));
        genotypeData = inputData.(genotypes{iGenotype}); 
        
        % Run analysis for high-firing cells only
        FRdata = genotypeData.highFiring;
        for iAnimal = 1:length(FRdata); 
            % Skip if empty
            if isempty(FRdata{iAnimal}) == 1; 
                continue
            else
                [~,n] = size(FRdata{iAnimal});
                for iCluster = 1:n;
                    % Skip if empty
                    if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                        display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                        continue
                    else
                        display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                        
                        % Assign variables based on running direction
                        directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original);
                        for iDir = 1:length(directions);
                            % Extract data based on running direction
                            variables = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir), 'spatialMetrics');
                            rateMap = variables.fullMap; timeMap = variables.timeMap; 
                            [m, n] = size(rateMap); 
                            
                            % Loop through trials
                            PFnumber = NaN(1,m); info = NaN(1,m); spars = NaN(1,m); 
                            for iTrial = 1:m;
                                % Get the size and number of the fields
                                [~, tempPFsize, PFnumber(iTrial)] = getPlaceFields_v1_20250425(rateMap(iTrial,:), thresholds);
                                for iField = 1:length(tempPFsize); 
                                    PFsize{iField}(iTrial) = tempPFsize(iField); 
                                end
                                % Get the position PDF for that trial
                                posPDF = timeMap(iTrial,:)/nansum(nansum(timeMap(iTrial,:))); 
                                % Get the information and sparsity
                                [info(iTrial), spars(iTrial), ~] = getMapStats(rateMap(iTrial,:), posPDF);
                            end
                            
                            % Save the output
                            if strcmp(directions(iDir), 'cw') == 1; 
                                outputData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).spatialMetrics.controls.byTrial.PFsize.cw = PFsize;
                                outputData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).spatialMetrics.controls.byTrial.PFnumber.cw = PFnumber;
                                outputData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).spatialMetrics.controls.byTrial.info.cw = info;
                                outputData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).spatialMetrics.controls.byTrial.sparsity.cw = spars;
                                fieldnames(outputData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).spatialMetrics)
                            end
                            % Save the output
                            if strcmp(directions(iDir), 'ccw') == 1; 
                                outputData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).spatialMetrics.controls.byTrial.PFsize.ccw = PFsize;
                                outputData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).spatialMetrics.controls.byTrial.PFnumber.ccw = PFnumber;
                                outputData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).spatialMetrics.controls.byTrial.info.ccw = info;
                                outputData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).spatialMetrics.controls.byTrial.sparsity.ccw = spars;
                            end
                        end
                    end
                end
            end
        end
    end
end
