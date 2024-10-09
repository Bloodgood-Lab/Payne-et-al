function [data, settings] = getPhasePrecession_v2_20241009(data, settings, processedDataPath)
    % Tries multiple methods for phase precession to find best approach
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
    
    data.cellData = data; data = rmfield(data, 'WT'); data = rmfield(data, 'KO');
    %% Step 1: Get the phase precession
    for iGenotype = 1:length(fieldnames(data.cellData));
        genotypes = fieldnames(data.cellData); 
        genotypeData = data.cellData.(genotypes{iGenotype}); 
        populationSlopes = []; rSquaredPopulation = []; populationCorrelation = [];
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
                        directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode);
                        for iDir = 1:length(directions);
                            % Extract variables based on running direction
                            outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir), 'phasePrecession');
                            spkPhs = outputData.spkPhs; spkPos = outputData.spkPos; binnedSpkPos = outputData.binnedSpkPos;
                            spikeTimes = outputData.spikesByDirection;
                            
                            % Loop through fields
                            numFieldsToAnalyze = whichField(settings.phasePrecession.fieldsToAnalyze, spkPhs);
                            slopeMedian = []; rSquared = []; rSquaredAllTrials = {}; slope = {}; spkPhsInput = {}; spkPosInput = {};
                            for iField = 1:numFieldsToAnalyze;
                                % Loop through all the trials
                                for iTrial = 1:length(spkPhs{iField});
                                    % If there are enough spatial bins
                                    if nanmax(binnedSpkPos{iField}{iTrial})-nanmin(binnedSpkPos{iField}{iTrial}) < settings.phasePrecession.spatialBinThreshold; 
                                        slope{iField}(iTrial) = NaN; 
                                        rSquaredAllTrials{iField}(iTrial) = NaN;
                                        continue; 
                                    else
                                        % If the number of spikes is more than threshold
                                        % and the time between spikes is within the threshold
                                        if length(spkPhs{iField}{iTrial}) >= settings.phasePrecession.spikeThreshold && ...
                                                max(abs(diff(spikeTimes{iField}{iTrial}))) <= settings.phasePrecession.ISIthreshold && ...
                                                max(spikeTimes{iField}{iTrial}) - min(spikeTimes{iField}{iTrial}) >= settings.phasePrecession.timeRange;
                                            
                                            % First approach: 
                                            % Get the phase precession using a circular slope-fitting method
                                            % described in Schmidt, 2009
                                            
                                            
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
                                            
                                            if strcmp(settings.phasePrecession.circularity, 'shift') == 1;
                                                % To account for circularity in the data, find the phase
                                                % shift that accounts for the best correlation between
                                                % position and spike theta phases
                                                spkPhsDegrees = rad2deg(spkPhs{iField}{iTrial});
                                                phasesToShiftBy = [0:180/36:360];
                                                correlation = [];
                                                for iShift = 1:length(phasesToShiftBy);
                                                    testShiftedPhs = circshift(spkPhsDegrees, phasesToShiftBy(iShift));
                                                    R = corrcoef(spkPosInput{iField}{iTrial}, testShiftedPhs); correlation(iShift) = R(2);
                                                end
                                                [~, maxInd] = min(correlation);
                                                spkPhsInput{iField}{iTrial} = circshift(spkPhsDegrees, phasesToShiftBy(maxInd));
                                                spkPhsInput{iField}{iTrial} = deg2rad(spkPhsInput{iField}{iTrial}) + pi;
                                            elseif strcmp(settings.phasePrecession.circularity, 'double') == 1;
                                                spkPosInput{iField}{iTrial} = [spkPosInput{iField}{iTrial}; spkPosInput{iField}{iTrial}];
                                                spkPhsInput{iField}{iTrial} = [spkPhs{iField}{iTrial}+pi; spkPhs{iField}{iTrial}+3*pi];
                                            elseif strcmp(settings.phasePrecession.circularity, 'none') == 1;
                                                spkPosInput{iField}{iTrial} = spkPosInput{iField}{iTrial};
                                                spkPhsInput{iField}{iTrial} = spkPhs{iField}{iTrial};
                                                R = corrcoef(spkPosInput{iField}{iTrial}, spkPhs{iField}{iTrial}); correlation = R(2); maxInd = 1; 
                                            end
                                            
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            %%%%%%% DEBUGGING STEPS %%%%%%%
                                            %figure(1); subplot(2,1,2); 
                                            %plot(spkPosInput{iField}{iTrial}, 'o-k'); 
                                            %xlim([0, length(spkPos{iField}{iTrial})+1]);
                                            %pause;
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            
                                            % Get the phase precession for that trial
                                            [cir, lin] = thetaPrecess(spkPhsInput{iField}{iTrial}, spkPosInput{iField}{iTrial}, settings.phasePrecession.slopeRange); 
                                            if strcmp(settings.phasePrecession.fit, 'circular') == 1;
                                                %y1 = [cir.Phi0, cir.Phi0+cir.Alpha];
                                                trialSlope = cir.Alpha; 
                                                trialOffset = cir.Phi0;
                                                r2 = (cir.Coeff)^2;
                                            elseif strcmp(settings.phasePrecession.fit, 'linear') == 1;
                                                %y1 = [lin.Phi0, lin.Phi0+lin.Alpha];
                                                trialSlope = lin.Alpha;
                                                trialOffset = lin.Phi0;
                                                r2 = (lin.r)^2;
                                            end
                                            if isempty(cir) == 0;
                                                fitInfo{iField}(iTrial).cir = cir; fitInfo{iField}(iTrial).lin = lin;
                                            else 
                                                fitInfo{iField}(iTrial).cir = NaN; fitInfo{iField}(iTrial).lin = NaN;
                                            end

                                            % If the slope is below the significance
                                            % threshold, save the data
                                            if cir.pValue < settings.phasePrecession.significanceThreshold; 
                                                %trialSlope = (y1(2) - y1(1))/(max(spkPosInput{iField}{iTrial})-min(spkPosInput{iField}{iTrial}));
                                                slope{iField}(iTrial) = trialSlope;
                                                offset{iField}(iTrial) = trialOffset; 
                                                rSquaredAllTrials{iField}(iTrial) = r2;
                                                if strcmp(settings.phasePrecession.circularity, 'shift') == 1 || ...
                                                        strcmp(settings.phasePrecession.circularity, 'none') == 1
                                                    phasePosCorrelation{iField}(iTrial) = correlation(maxInd);
                                                else
                                                    phasePosCorrelation{iField}(iTrial) = NaN;
                                                end
                                            else
                                                slope{iField}(iTrial) = NaN;
                                                offset{iField}(iTrial) = NaN; 
                                                rSquaredAllTrials{iField}(iTrial) = NaN;
                                                phasePosCorrelation{iField}(iTrial) = NaN;
                                            end
                                            
                                        else
                                            spkPhsInput{iField}{iTrial} = []; 
                                            spkPosInput{iField}{iTrial} = []; 
                                            slope{iField}(iTrial) = NaN;
                                            offset{iField}(iTrial) = NaN; 
                                            rSquaredAllTrials{iField}(iTrial) = NaN;
                                            fitInfo{iField}(iTrial).cir = NaN; 
                                            fitInfo{iField}(iTrial).lin = NaN;
                                            phasePosCorrelation{iField}(iTrial) = NaN;
                                        end
                                    end
                                    
                                end
                                
                                % If there are enough trials with slopes
                                % calculated, get the median of all slopes
                                if sum(~isnan(slope{iField})) >= settings.phasePrecession.trialThreshold; 
                                    slopeMedian(iField) = nanmedian(slope{iField});
                                    rSquared(iField) = nanmean(rSquaredAllTrials{iField});
                                else 
                                    slopeMedian(iField) = NaN; 
                                    rSquared(iField) = NaN;
                                end
                                
                                % Assign output data
                                if strcmp(directions(iDir), 'cw') == 1;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.phsInput.cw = spkPhsInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.posInput.cw = spkPosInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.allSlopes.cw = slope;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSlope.cw = slopeMedian;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.fitInfo.cw = fitInfo;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.rSquared.cw = rSquaredAllTrials;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.correlation.cw = phasePosCorrelation;
                                elseif strcmp(directions(iDir), 'ccw') == 1; 
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.phsInput.ccw = spkPhsInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.posInput.ccw = spkPosInput;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.allSlopes.ccw = slope;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.medianSlope.ccw = slopeMedian;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.fitInfo.ccw = fitInfo;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.rSquared.ccw = rSquaredAllTrials;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).phasePrecession.correlation.ccw = phasePosCorrelation;
                                end
                                allTrialsSlopes = cell2mat(slope); 
                                allTrialsOffsets = cell2mat(offset);
                                allTrialsCorrelation = cell2mat(phasePosCorrelation);
                            end
                            
                            % Organize the population data and analysis
                            populationSlopes = [populationSlopes, slopeMedian];
                            rSquaredPopulation = [rSquaredPopulation, rSquared]; 
                            populationAllTrialSlopes = [populationAllTrialSlopes, allTrialsSlopes]; 
                            populationAllTrialOffsets = [populationAllTrialOffsets, allTrialsOffsets]; 
                            populationCorrelation = [populationCorrelation, allTrialsCorrelation];
                        end
                    end
                end
            end
        end
        data.populationData(iGenotype).phasePrecessionSlopes = populationSlopes;
        data.populationData(iGenotype).phasePrecessionRSquared = rSquaredPopulation;
        data.populationData(iGenotype).phasePrecessionAllTrialSlopes = populationAllTrialSlopes;
        data.populationData(iGenotype).phasePrecessionAllTrialOffsets = populationAllTrialOffsets;
        data.populationData(iGenotype).phasePositionCorrelation = populationCorrelation;
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


    