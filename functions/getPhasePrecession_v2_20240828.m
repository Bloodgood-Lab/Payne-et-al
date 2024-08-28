function data = getPhasePrecession_v2_20240828(data, settings, processedDataPath)
    % Calculates a linear fit to get the phase precession. This approach
    % was discarded in favor of the approach described in Robert Schmidt, 
    % 2009 (see 'getPhasePrecession_v1_20240806) which better accounts for
    % the circularity of the data
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
    
    if settings.phasePrecession.spikeThreshold < 2; 
        error('Spike minimum must be at least 2, please change settings');
    end
    
    data.cellData = data; data = rmfield(data, 'WT'); data = rmfield(data, 'KO');
    %% Step 1: Get the phase precession
    for iGenotype = 1%:length(fieldnames(data.cellData));
        genotypes = fieldnames(data.cellData); 
        genotypeData = data.cellData.(genotypes{iGenotype}); 
        populationSlopes = []; 
        
        % Run analysis for high-firing cells only
        FRdata = genotypeData.highFiring;
        for iAnimal = 1%1:length(FRdata); 
            % Skip if empty
            if isempty(FRdata{iAnimal}) == 1; 
                continue
            else
                [~,n] = size(FRdata{iAnimal});
                for iCluster = 2%1:n;
                    % Skip if empty
                    if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                        display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                        continue
                    else
                        display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                        
                        % Assign variables based on running direction
                        directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode);
                        for iDir = 1:length(directions);
                            if strcmp(directions(iDir), 'cw') == 1; 
                                spkPhs = FRdata{iAnimal}(iCluster).theta.phases.cw; 
                                spkPos = FRdata{iAnimal}(iCluster).inField.inFieldSpkPos.cw; 
                                binnedSpkPos = FRdata{iAnimal}(iCluster).inField.inFieldBinnedSpkPos.cw; 
                            elseif strcmp(directions(iDir), 'ccw') == 1; 
                                spkPhs = FRdata{iAnimal}(iCluster).theta.phases.ccw; 
                                spkPos = FRdata{iAnimal}(iCluster).inField.inFieldSpkPos.ccw; 
                                binnedSpkPos = FRdata{iAnimal}(iCluster).inField.inFieldBinnedSpkPos.ccw; 
                            end
                        
                            % Loop through fields
                            if strcmp(settings.phasePrecession.fieldsToAnalyze, 'all fields') == 1;
                                numFieldsToAnalyze = length(spkPhs); 
                            elseif strcmp(settings.phasePrecession.fieldsToAnalyze, 'best field') == 1;
                                numFieldsToAnalyze = 1; 
                            end
                            slopeMedian = []; subplotCount = 1; 
                            for iField = 1:numFieldsToAnalyze;
                                % Loop through all the trials
                                slope{iField} = NaN(1,length(spkPhs{iField}));
                                spkPosInput{iField} = cell(1, length(spkPos{iField}));
                                spkPhsInput{iField} = cell(1, length(spkPhs{iField}));
                                y_fit{iField} = cell(1,length(spkPhs{iField})); 
                                for iTrial = 1:length(spkPhs{iField});
                                    % If there were spikes in-field that trial
                                    if length(spkPhs{iField}{iTrial}) >= settings.phasePrecession.spikeThreshold;
                                        % Convert phases to degrees for
                                        % ease in working with the data
                                        spkPhsDegrees = rad2deg(spkPhs{iField}{iTrial});
                                    
                                        % Determine whether to use
                                        % binned or unbinned position
                                        if strcmp(settings.phasePrecession.positionType, 'unbinned') == 1
                                            spkPosToUse = spkPos{iField}{iTrial};
                                            % If it's a CW trial, flip it
                                            % so that the position is in
                                            % the direction of running
                                            if strcmp(directions(iDir), 'cw') == 1;
                                                spkPosToUse = fliplr(spkPos{iField}{iTrial});
                                            end
                                        elseif strcmp(settings.phasePrecession.positionType, 'binned') == 1
                                            spkPosToUse = binnedSpkPos{iField}{iTrial};
                                            % If it's a CW trial, flip it
                                            % so that the position is in
                                            % the direction of running
                                            if strcmp(directions(iDir), 'cw') == 1;
                                                spkPosToUse = fliplr(spkPos{iField}{iTrial});
                                            end
                                        end
                                            
                                        % Normalize position between 0
                                        % and 1
                                        spkPosInput{iField}{iTrial} = (spkPosToUse-nanmin(spkPosToUse))/(nanmax(spkPosToUse)-nanmin(spkPosToUse));

                                        % To account for circularity in
                                        % the data, find the phase
                                        % shift that accounts for the
                                        % best correlation between
                                        % position and spike theta
                                        % phases
                                        phasesToShiftBy = [0:180/36:180];
                                        correlation = [];
                                        for iShift = 1:length(phasesToShiftBy);
                                            testShiftedPhs = circshift(spkPhsDegrees, phasesToShiftBy(iShift));
                                            R = corrcoef(spkPosInput{iField}{iTrial}, testShiftedPhs); correlation(iShift) = R(2);
                                        end
                                        [~, maxInd] = max(abs(correlation));
                                        spkPhsInput{iField}{iTrial} = circshift(spkPhsDegrees, phasesToShiftBy(maxInd));
                                        
                                        % If there are enough spatial bins
                                        if nanmax(binnedSpkPos{iField}{iTrial})-nanmin(binnedSpkPos{iField}{iTrial}) < settings.phasePrecession.spatialBinThreshold;                                             slope{iField}(iTrial) = NaN;
                                            continue;
                                        else
                                            % Get the phase precession for that trial
                                            p{iField}{iTrial} = polyfit(spkPosInput{iField}{iTrial}, spkPhsInput{iField}{iTrial}, 1);  % p(1) is the slope, p(2) is the intercept
                                            y_fit{iField}{iTrial} = polyval(p{iField}{iTrial}, spkPosInput{iField}{iTrial});
                                            slope{iField}(iTrial) = p{iField}{iTrial}(1);
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
                                if strcmp(directions(iDir), 'cw') == 1;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.allSlopes.cw = slope;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSlope.cw = slopeMedian;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.posInput.cw = spkPosInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.phsInput.cw = spkPhsInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.yFit.cw = y_fit;
                                elseif strcmp(directions(iDir), 'ccw') == 1; 
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.allSlopes.ccw = slope;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSlope.ccw = slopeMedian;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.phsInput.ccw = spkPosInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.phsInput.ccw = spkPhsInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.yFit.ccw = y_fit;
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
    saveFile_v1_20240718(processedDataPath, data, settings, 'phasePrecession')
    
    
    