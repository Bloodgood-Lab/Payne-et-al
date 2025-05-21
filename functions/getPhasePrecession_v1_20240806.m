function [data, settings] = getPhasePrecession_v1_20240806(data, settings, processedDataPath)
    % Relies on method described in "Robert Schmidt, 2009, Single-Trial 
    % Place Precession in the Hippocampus" to get the phase precession
    % Written by Anja Payne
    % Last Modified: 10/09/2024

    % Inputs:
    %   1) data: the matlab structure where the in-field theta phases
    %      by-trial are saved
    %   2) settings: the settings file where the settings for the phase
    %      precession calculations are saved
    
    % Outputs:
    %   1) data: the structure with the phase precession slopes appended

    % Steps:
    %   1) Loop through genotypes, animals, and cells and get the phase
    %      precession
    %   2) Ask the user if they want to save the newly generated data
    %   3) Run the statistics on the data and output the p-value
    
    if isfield(data, 'WT')
        data.cellData = data; data = rmfield(data, 'WT'); data = rmfield(data, 'KO');
    end
    %% Step 1: Get the phase precession
    for iGenotype = 1:length(fieldnames(data.cellData));
        genotypes = fieldnames(data.cellData); 
        genotypeData = data.cellData.(genotypes{iGenotype}); 
        populationSlopes = []; populationMedianSizes = []; populationMeanSizes = [];
        population_thresh_spatialBin = []; population_thresh_spikeNum = []; 
        population_thresh_ISI = []; population_thresh_timeRange = []; 
        population_thresh_p = []; population_thresh_trialNum = []; 
        rSquaredPopulation = [];
        populationAllTrialSlopes = []; populationAllTrialOffsets = []; 
        
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
                            % Extract variables based on running direction
                            outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir), 'phasePrecession');
                            spkPhs = outputData.spkPhs; spkPos = outputData.spkPos; binnedSpkPos = outputData.binnedSpkPos;
                            spikeTimes = outputData.spikesByDirection;
                            
                            % Loop through fields
                            numFieldsToAnalyze = whichField(settings.phasePrecession.fieldsToAnalyze, spkPhs);
                            slopeMedian = []; sizeMedian = []; rSquared = []; rSquaredAllTrials = {}; 
                            slope = {}; spkPhsInput = {}; spkPosInput = {}; fieldSize = {}; 
                            trialSlope = {}; trialOffset = {}; trialR2 = {}; trialPvalue = {}; 
                            for iField = 1:numFieldsToAnalyze;
                                % Loop through all the trials
                                allTrialsPosition = []; allTrialsPhases = []; fieldSize{iField} = []; 
                                for iTrial = 1:length(spkPhs{iField});
                                    % If there are enough spatial bins
                                    if nanmax(binnedSpkPos{iField}{iTrial})-nanmin(binnedSpkPos{iField}{iTrial}) < settings.phasePrecession.spatialBinThreshold; 
                                        slope{iField}(iTrial) = NaN; 
                                        rSquaredAllTrials{iField}(iTrial) = NaN;
                                        thresh_spatialBin = 0; 
                                        continue; 
                                    else
                                        thresh_spatialBin = 1; 
                                        % If the number of spikes is more than threshold
                                        % and the time between spikes is within the threshold
                                        if length(spkPhs{iField}{iTrial}) >= settings.phasePrecession.spikeThreshold && ...
                                                max(abs(diff(spikeTimes{iField}{iTrial}))) <= settings.phasePrecession.ISIthreshold && ...
                                                max(spikeTimes{iField}{iTrial}) - min(spikeTimes{iField}{iTrial}) >= settings.phasePrecession.timeRange;
                                            
                                            % Keep track of thresholds passed
                                            thresh_spikeNum = 1; 
                                            thresh_ISI = 1; 
                                            thresh_timeRange = 1; 
                                            
                                            % Get the position information
                                            if strcmp(settings.phasePrecession.positionType, 'unbinned') == 1
                                                spkPosInput{iField}{iTrial} = spkPos{iField}{iTrial}-nanmin(spkPos{iField}{iTrial});
                                            elseif strcmp(settings.phasePrecession.positionType, 'binned') == 1
                                                spkPosInput{iField}{iTrial} = binnedSpkPos{iField}{iTrial}-min(binnedSpkPos{iField}{iTrial});
                                            end
                                            % For cw trials, values will be decreasing since scale is
                                            % absolute. To fit slope accurately, they need to be increasing. 
                                            if strcmp(directions(iDir), 'cw') == 1; 
                                                spkPosInput{iField}{iTrial} = abs(spkPosInput{iField}{iTrial}-nanmax(spkPosInput{iField}{iTrial})); 
                                            end
                                            if strcmp(settings.phasePrecession.normalized, 'yes') == 1; 
                                                spkPosInput{iField}{iTrial} = spkPosInput{iField}{iTrial}/max(spkPosInput{iField}{iTrial});
                                            end
                                            
                                            % Calculate the phase precession depending on the user settings
                                            inputData.spkPhs = spkPhs{iField}{iTrial}; 
                                            inputData.spkPos = spkPosInput{iField}{iTrial};
                                            
                                            outputData = calculatePhasePrecession_v1_20241020(inputData, settings); 
                                            
                                            spkPhsInput{iField}{iTrial} = outputData.spkPhs; 
                                            spkPosInput{iField}{iTrial} = outputData.spkPos; 
                                            fieldSize{iField}(iTrial) = nanmax(outputData.spkPos) - nanmin(outputData.spkPos);
                                            trialSlope{iField}(iTrial) = outputData.trialSlope; 
                                            trialOffset{iField}(iTrial) = outputData.trialOffset; 
                                            trialR2{iField}(iTrial) = outputData.r2; 
                                            trialPvalue{iField}(iTrial) = outputData.p; 
                                            
                                            % If the slope is below the significance threshold, save the data
                                            if trialPvalue{iField}(iTrial) > settings.phasePrecession.significanceThreshold; 
                                                fieldSize{iField}(iTrial) = NaN; 
                                                trialSlope{iField}(iTrial) = NaN;
                                                trialOffset{iField}(iTrial) = NaN; 
                                                trialR2{iField}(iTrial) = NaN;
                                                trialPvalue{iField}(iTrial) = NaN; 
                                                
                                                % Keep track of thresholds passed
                                                thresh_p = 0;
                                                
                                            else
                                                thresh_p = 1; 
                                            end
                                            
                                        else
                                            spkPhsInput{iField}{iTrial} = []; 
                                            spkPosInput{iField}{iTrial} = []; 
                                            trialSlope{iField}(iTrial) = NaN;
                                            trialOffset{iField}(iTrial) = NaN; 
                                            trialR2{iField}(iTrial) = NaN;
                                            trialPvalue{iField}(iTrial) = NaN;
                                            fieldSize{iField}(iTrial) = NaN; 
                                            
                                            % Keep track of thresholds passed
                                            thresh_spikeNum = 1; 
                                            thresh_ISI = 1; 
                                            thresh_timeRange = 1; 
                                        end
                                    end
                                end
                                
                                % If there are enough trials with slopes calculated, get the median of all slopes
                                if sum(~isnan(trialSlope{iField})) >= settings.phasePrecession.trialThreshold; 
                                    slopeMedian(iField) = nanmedian(trialSlope{iField});
                                    [~, idx] = nanmin(abs(trialSlope{iField} - slopeMedian(iField)));
                                    sizeMedian(iField) = nanmedian(fieldSize{iField})
                                    sizeMean(iField) = nanmean(fieldSize{iField});
                                    rSquared(iField) = nanmean(trialR2{iField});
                                    
                                    % Keep track of thresholds passed
                                    thresh_trialNum = 1; 
                                else 
                                    slopeMedian(iField) = NaN; 
                                    rSquared(iField) = NaN;
                                    sizeMedian(iField) = NaN; 
                                    sizeMean(iField) = NaN; 
                                    
                                    % Keep track of thresholds passed
                                    thresh_trialNum = 0; 
                                end
                                
                                % Assign output data
                                if strcmp(directions(iDir), 'cw') == 1;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.phsInput.cw = spkPhsInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.posInput.cw = spkPosInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.allSlopes.cw = trialSlope;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.offsets.cw = trialOffset;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.rSquared.cw = trialR2;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.pValues.cw = trialPvalue;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSlope.cw = slopeMedian;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSize.cw = sizeMedian;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.meanR2.cw = rSquared;
                                elseif strcmp(directions(iDir), 'ccw') == 1; 
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.phsInput.ccw = spkPhsInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.posInput.ccw = spkPosInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.allSlopes.ccw = trialSlope;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.offsets.ccw = trialOffset;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.rSquared.ccw = trialR2;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.pValues.ccw = trialPvalue;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSlope.ccw = slopeMedian;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSize.ccw = sizeMedian;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.meanR2.ccw = rSquared;
                                end
                                
                                allTrialsSlopes = cell2mat(trialSlope); 
                                allTrialsOffsets = cell2mat(trialOffset);

                            end
                            
                            % Organize the population data and analysis
                            populationSlopes = [populationSlopes, slopeMedian];
                            rSquaredPopulation = [rSquaredPopulation, rSquared]; 
                            populationAllTrialSlopes = [populationAllTrialSlopes, allTrialsSlopes]; 
                            populationAllTrialOffsets = [populationAllTrialOffsets, allTrialsOffsets]; 
                            populationMedianSizes = [populationMedianSizes, sizeMedian]; 
                            populationMeanSizes = [populationMeanSizes, sizeMean]; 
                            population_thresh_spatialBin = [population_thresh_spatialBin, thresh_spatialBin]; 
                            population_thresh_spikeNum = [population_thresh_spikeNum, thresh_spikeNum]; 
                            population_thresh_ISI = [population_thresh_ISI, thresh_ISI];
                            population_thresh_timeRange = [population_thresh_timeRange, thresh_timeRange]; 
                            population_thresh_p = [population_thresh_p, thresh_p];
                            population_thresh_trialNum = [population_thresh_trialNum, thresh_trialNum]; 
                        end
                    end
                end
            end
        end
        data.populationData(iGenotype).phasePrecession.Slopes = populationSlopes;
        data.populationData(iGenotype).phasePrecession.RSquared = rSquaredPopulation;
        data.populationData(iGenotype).phasePrecession.AllTrialSlopes = populationAllTrialSlopes;
        data.populationData(iGenotype).phasePrecession.AllTrialOffsets = populationAllTrialOffsets;
        data.populationData(iGenotype).phasePrecession.MedianFieldSizes = populationMedianSizes;
        data.populationData(iGenotype).phasePrecession.MeanFieldSizes = populationMeanSizes;
        data.populationData(iGenotype).phasePrecession.thresholds.spatialBin = population_thresh_spatialBin;
        data.populationData(iGenotype).phasePrecession.thresholds.ISI = population_thresh_ISI;
        data.populationData(iGenotype).phasePrecession.thresholds.timeRange = population_thresh_timeRange;
        data.populationData(iGenotype).phasePrecession.thresholds.pValue = population_thresh_p;
        data.populationData(iGenotype).phasePrecession.thresholds.trialNum = population_thresh_trialNum;
        
    end
    
    %% Step 2: Save
    settings = saveFile_v1_20240718(processedDataPath, data, settings, 'phasePrecession')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function numField = whichField(fieldsToAnalyze, spikes)  
        % Loop through fields
        if strcmp(fieldsToAnalyze, 'all fields') == 1;
            numField = length(spikes); 
        elseif strcmp(fieldsToAnalyze, 'best field') == 1;
            numField = 1; 
        end
    end


    