function [data, settings] = getPhasePrecession_v1_20240806(data, settings, processedDataPath)
    % Relies on method described in "Robert Schmidt, 2009, Single-Trial 
    % Place Precession in the Hippocampus" to get the phase precession
    % Written by Anja Payne
    % Last Modified: 08/28/2024

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
    
    data.cellData = data; data = rmfield(data, 'WT'); data = rmfield(data, 'KO');
    %% Step 1: Get the phase precession
    for iGenotype = 1%:length(fieldnames(data.cellData));
        genotypes = fieldnames(data.cellData); 
        genotypeData = data.cellData.(genotypes{iGenotype}); 
        populationSlopes = []; 
        
        % Run analysis for high-firing cells only
        FRdata = genotypeData.highFiring;
        for iAnimal = 1%:length(FRdata); 
            % Skip if empty
            if isempty(FRdata{iAnimal}) == 1; 
                continue
            else
                [~,n] = size(FRdata{iAnimal});
                for iCluster = 4%1:n;
                    % Skip if empty
                    if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                        display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                        continue
                    else
                        display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                        
                        % Assign variables based on running direction
                        directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode);
                        for iDir = 1:length(directions);
                            % Extract variables based on running direction
                            outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir), 'phasePrecession');
                            spkPhs = outputData.spkPhs; spkPos = outputData.spkPos; binnedSpkPos = outputData.binnedSpkPos;
                            spikeTimes = outputData.spikesByDirection;
                            
                            
                            % Loop through fields
                            numFieldsToAnalyze = whichField(settings.phasePrecession.fieldsToAnalyze, spkPhs);
                            slopeMedian = []; slope = {}; spkPhsInput = {}; spkPosInput = {};
                            for iField = 1:numFieldsToAnalyze;
                                % Loop through all the trials
                                for iTrial = 1:length(spkPhs{iField});
                                    % If there are enough spatial bins
                                    if nanmax(binnedSpkPos{iField}{iTrial})-nanmin(binnedSpkPos{iField}{iTrial}) < settings.phasePrecession.spatialBinThreshold; 
                                        slope{iField}(iTrial) = NaN; 
                                        continue; 
                                    else
                                        % If the number of spikes is more than threshold
                                        % and the time between spikes is within the threshold
                                        if length(spkPhs{iField}{iTrial}) >= settings.phasePrecession.spikeThreshold && ...
                                                max(abs(diff(spikeTimes{iField}{iTrial}))) <= settings.phasePrecession.ISIthreshold; 
                                            % Get the phase precession for that trial
                                            spkPhsInput{iField}{iTrial} = [spkPhs{iField}{iTrial}+pi; spkPhs{iField}{iTrial}+3*pi]; 
                                            if strcmp(settings.phasePrecession.positionType, 'unbinned') == 1
                                                spkPosInput{iField}{iTrial} = [spkPos{iField}{iTrial}-min(spkPos{iField}{iTrial}); spkPos{iField}{iTrial}-min(spkPos{iField}{iTrial})];
                                            elseif strcmp(settings.phasePrecession.positionType, 'binned') == 1
                                                spkPosInput{iField}{iTrial} = [binnedSpkPos{iField}{iTrial}-min(binnedSpkPos{iField}{iTrial}); binnedSpkPos{iField}{iTrial}-min(binnedSpkPos{iField}{iTrial})];
                                            end
                                            % For cw trials, values will be decreasing since scale is
                                            % absolute. To fit slope accurately, they need to be increasing. 
                                            if strcmp(directions(iDir), 'cw') == 1; 
                                                spkPosInput{iField}{iTrial} = abs(spkPosInput{iField}{iTrial}-nanmax(spkPosInput{iField}{iTrial})); 
                                            end;
                                            spkPosInput{iField}{iTrial} = spkPosInput{iField}{iTrial}/max(spkPosInput{iField}{iTrial}); 
                                            [cir, lin] = thetaPrecess(spkPhsInput{iField}{iTrial}, spkPosInput{iField}{iTrial}, settings.phasePrecession.slopeRange); 
                                            y1 = [cir.Phi0, cir.Phi0+cir.Alpha];
                                            if isempty(cir) == 0;
                                                fitInfo{iField}(iTrial).cir = cir; fitInfo{iField}(iTrial).lin = lin;
                                            else 
                                                fitInfo{iField}(iTrial).cir = NaN; fitInfo{iField}(iTrial).lin = NaN;
                                            end

                                            % If the slope is below the significance
                                            % threshold, save the data
                                            if cir.pValue < settings.phasePrecession.significanceThreshold; 
                                                trialSlope = (y1(2) - y1(1))/(max(spkPos{iField}{iTrial})-min(spkPos{iField}{iTrial}));
                                                slope{iField}(iTrial) = trialSlope;
                                            else
                                                slope{iField}(iTrial) = NaN; 
                                            end
                                            
                                        else
                                            spkPhsInput{iField}{iTrial} = []; 
                                            spkPosInput{iField}{iTrial} = []; 
                                            slope{iField}(iTrial) = NaN;
                                        end
                                    end
                                end
                                
                                % If there are enough trials with slopes
                                % calculated, get the median of all slopes
                                if sum(~isnan(slope{iField})) >= settings.phasePrecession.trialThreshold; 
                                    slopeMedian(iField) = nanmedian(slope{iField});
                                else 
                                    slopeMedian(iField) = NaN; 
                                end
                                
                                % Assign output data
                                if strcmp(directions(iDir), 'cw') == 1;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.phsInput.cw = spkPhsInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.posInput.cw = spkPosInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.allSlopes.cw = slope;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSlope.cw = slopeMedian;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.fitInfo.cw = fitInfo;
                                elseif strcmp(directions(iDir), 'ccw') == 1; 
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.phsInput.ccw = spkPhsInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.posInput.ccw = spkPosInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.allSlopes.ccw = slope;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSlope.ccw = slopeMedian;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.fitInfo.ccw = fitInfo;
                                end
                            end
                            
                            % Organize the population data and analysis
                            populationSlopes = [populationSlopes, slopeMedian];
                        end
                    end
                end
            end
        end
        data.populationData(iGenotype).phasePrecessionSlopes = populationSlopes;
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


    