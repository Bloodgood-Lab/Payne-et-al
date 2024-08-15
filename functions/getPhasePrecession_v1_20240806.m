function data = getPhasePrecession_v1_20240806(data, settings, processedDataPath)
    % Place Precession in the Hippocampus))
    % Gets the theta angle of spikes
    % Written by Anja Payne
    % Last Modified: 08/15/2024

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
    
    %% Step 1: Get the phase precession
    for iGenotype = 1:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype}); 
        
        % Run analysis for high-firing cells only
        FRdata = genotypeData.highFiring;
        for iAnimal = 1:3%:length(FRdata); 
            % Skip if empty
            if isempty(FRdata{iAnimal}) == 1; 
                continue
            else
                [~,n] = size(FRdata{iAnimal});
                for iCluster = 1%:n;
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
                            elseif strcmp(directions(iDir), 'ccw') == 1; 
                                spkPhs = FRdata{iAnimal}(iCluster).theta.phases.ccw; 
                                spkPos = FRdata{iAnimal}(iCluster).inField.inFieldSpkPos.ccw; 
                            end
                        end
                        
                        % Loop through all the fields
                        slopeMedian = [];
                        for iField = 1:length(spkPhs);
                            % Loop through all the trials
                            slope{iField} = NaN(1,length(spkPhs{iField}));
                            for iTrial = 1:length(spkPhs{iField});
                                % If there are enough spatial bins
                                if nanmax(spkPos{iField}{iTrial})-nanmin(spkPos{iField}{iTrial}) < settings.phasePrecession.spatialBinThreshold; 
                                    continue; 
                                else
                                    % If there were spikes in-field that trial
                                    if isempty(spkPhs{iField}{iTrial}) == 0; 

                                        % Get the phase precession for that trial
                                        spkPhsInput = [spkPhs{iField}{iTrial}+pi; spkPhs{iField}{iTrial}+3*pi]; 
                                        spkPosInput = [spkPos{iField}{iTrial}; spkPos{iField}{iTrial}];
                                        [cir, lin] = thetaPrecess(spkPhsInput, spkPosInput-min(spkPosInput), settings.phasePrecession.slopeRange); 
                                        y1 = [cir.Phi0, cir.Phi0+cir.Alpha];

                                        % If the slope is below the significance
                                        % threshold, save the data
                                        if cir.pValue < settings.phasePrecession.significanceThreshold; 
                                            trialSlope = (y1(2) - y1(1))/(max(spkPos{iField}{iTrial})-min(spkPos{iField}{iTrial}));
                                            slope{iField}(iTrial) = trialSlope;
                                        else
                                            slope{iField}(iTrial) = NaN; 
                                        end
                                    else
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

                            if strcmp(directions(iDir), 'cw') == 1;
                                data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.allSlopes.cw = slope;
                                data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSlope.cw = slopeMedian;
                            elseif strcmp(directions(iDir), 'ccw') == 1; 
                                data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.allSlopes.ccw = slope;
                                data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSlope.ccw = slopeMedian;
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% Step 2: Save
    saveFile_v1_20240718(processedDataPath, data, settings, 'phasePrecession') 
    
    
    