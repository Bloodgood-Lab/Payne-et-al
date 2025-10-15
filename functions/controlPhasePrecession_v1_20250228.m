function data = controlPhasePrecession_v1_20250228(data, settings)
    % Generates plots related to phase precession analysis
    % Written by Anja Payne
    % Last Modified: 03/11/2025
    
    % Inputs:
    %   1) data: the matlab structure where the in-field theta phases
    %      by-trial are saved
    %   2) settings: the settings file where the settings for the phase
    %      precession calculations are saved    
    
    % Steps/Figures:
    %   1) Get the average size of the field for each cell that was
    %      included in the phase precession slope analysis
    %   2) Match the distributions of the place field sizes then compare 
    %      the slopes
    %   3) Match the distributions of the average place field sizes then
    %      compare the slopes
    %   4) Bootstrap the cells to determine whether we have enough cells in
    %      each group
    
    
    %% Step 1: Get the average place field size
    for iGenotype = 1:length(fieldnames(data.cellData));
        genotypes = fieldnames(data.cellData); 
        genotypeData = data.cellData.(genotypes{iGenotype}); 
        allFieldSizes = []; 
        
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
                            % Loop through fields
                            fieldsToAnalyze = settings.phasePrecession.fieldsToAnalyze; 
                            if strcmp(fieldsToAnalyze, 'all fields') == 1;
                                numField = length(FRdata{iAnimal}(iCluster).spatialMetrics.spatialMetrics.PFsize.(directions{iDir})); 
                            elseif strcmp(fieldsToAnalyze, 'best field') == 1;
                                numField = 1; 
                            end
                            
                            for iField = 1:numField;
                                singleFieldSize = FRdata{iAnimal}(iCluster).spatialMetrics.PFsize.(directions{iDir})(iField);  
                                allFieldSizes = [allFieldSizes, singleFieldSize];
                            end
                        end
                    end
                end
            end
        end
        data.populationData(iGenotype).phasePrecession.control.averageFieldSizes = allFieldSizes;
    end
    
    %% Step 2: Match the distributions of the median place field sizes then compare the slopes
    for iGenotype = 1:length(data.populationData)
        
        % Sort by size and apply the same sorting to the slope vector
        [sortedSizes{iGenotype}, sortedSize_idx] = sort(data.populationData(iGenotype).phasePrecession.MedianFieldSizes);
        sortedSlopes{iGenotype} = data.populationData(iGenotype).phasePrecession.Slopes(sortedSize_idx);
        
        % Bin the size vector up
        edgeValues = [0:4:max(sortedSizes{1})];
        [sizeHistogramValues(iGenotype, :), ~] = histcounts(sortedSizes{iGenotype}, edgeValues);
        
    end

    % Initialize variables
    numIterations = 1;
    maxN = length(sizeHistogramValues);
    edgeValues = [0:4:max(sortedSizes{1})];

    % Preallocate results
    p_slope = zeros(numIterations, 1);
    p_size = zeros(numIterations, 1);

    % Loop through random iterations
    for iRandomIteration = 1:numIterations
        iRandomIteration
        DSsize_WT = []; DSsize_KO = []; 
        DSslope_WT = []; DSslope_KO = []; 

        % Loop through each bin
        for iHistBin = 1:maxN-1
            % Find indices of values within the current bin
            WTindex = find(sortedSizes{1} >= edgeValues(iHistBin) & sortedSizes{1} < edgeValues(iHistBin+1));
            KOindex = find(sortedSizes{2} >= edgeValues(iHistBin) & sortedSizes{2} < edgeValues(iHistBin+1));

            % Downsample WT if there are more WT values than KO
            if length(WTindex) > length(KOindex)
                randomWTindex = randperm(length(WTindex), length(KOindex));
                DSsize_WT = [DSsize_WT, sortedSizes{1}(WTindex(randomWTindex))];
                DSslope_WT = [DSslope_WT, sortedSlopes{1}(WTindex(randomWTindex))];
                DSsize_KO = [DSsize_KO, sortedSizes{2}(KOindex)];
                DSslope_KO = [DSslope_KO, sortedSlopes{2}(KOindex)];
            % Downsample KO if there are more KO values than WT
            elseif length(KOindex) > length(WTindex)
                randomKOindex = randperm(length(KOindex), length(WTindex));
                DSsize_KO = [DSsize_KO, sortedSizes{2}(KOindex(randomKOindex))];
                DSslope_KO = [DSslope_KO, sortedSlopes{2}(KOindex(randomKOindex))];
                DSsize_WT = [DSsize_WT, sortedSizes{1}(WTindex)];
                DSslope_WT = [DSslope_WT, sortedSlopes{1}(WTindex)];
            % If the sizes are equal, just append all values
            else
                DSsize_WT = [DSsize_WT, sortedSizes{1}(WTindex)];
                DSslope_WT = [DSslope_WT, sortedSlopes{1}(WTindex)];
                DSsize_KO = [DSsize_KO, sortedSizes{2}(KOindex)];
                DSslope_KO = [DSslope_KO, sortedSlopes{2}(KOindex)];
            end
        end

        % Perform KS test for slopes and sizes
        [~, p_slope(iRandomIteration)] = kstest2(DSslope_WT, DSslope_KO);
        [~, p_size(iRandomIteration)] = kstest2(DSsize_WT, DSsize_KO);
        
    end
    data.populationData(1).phasePrecession.control.matchDistribution.byMedian.dsSizes = DSsize_WT;
    data.populationData(2).phasePrecession.control.matchDistribution.byMedian.dsSizes = DSsize_KO;
    data.populationData(1).phasePrecession.control.matchDistribution.byMedian.dsSlopes = DSslope_WT;
    data.populationData(2).phasePrecession.control.matchDistribution.byMedian.dsSlopes = DSslope_KO;
    data.populationData(1).phasePrecession.control.matchDistribution.byMedian.sizePvalues = p_size;
    data.populationData(1).phasePrecession.control.matchDistribution.byMedian.slopePvalues = p_slope;
    
    %% Step 3: Match the distributions of the average place field sizes then compare the slopes
    sizeHistogramValues = []; 
    for iGenotype = 1:length(data.populationData)
        
        % Sort by size and apply the same sorting to the slope vector
        [sortedSizes{iGenotype}, sortedSize_idx] = sort(data.populationData(iGenotype).phasePrecession.MeanFieldSizes);
        sortedSlopes{iGenotype} = data.populationData(iGenotype).phasePrecession.Slopes(sortedSize_idx);
        
        % Bin the size vector up
        edgeValues = [0:4:max(sortedSizes{1})];
        [sizeHistogramValues(iGenotype, :), ~] = histcounts(sortedSizes{iGenotype}, edgeValues);
        
    end
    
    % Set random seed for reproducibility
    rng('shuffle');

    % Initialize variables
    numIterations = 10;
    maxN = length(sizeHistogramValues);
    edgeValues = [0:4:max(sortedSizes{1})];

    % Preallocate results
    p_slope = zeros(numIterations, 1);
    p_size = zeros(numIterations, 1);

    % Loop through random iterations
    for iRandomIteration = 1:numIterations
        iRandomIteration
        DSsize_WT = []; DSsize_KO = []; 
        DSslope_WT = []; DSslope_KO = []; 

        % Loop through each bin
        for iHistBin = 1:maxN-1
            % Find indices of values within the current bin
            WTindex = find(sortedSizes{1} >= edgeValues(iHistBin) & sortedSizes{1} < edgeValues(iHistBin+1));
            KOindex = find(sortedSizes{2} >= edgeValues(iHistBin) & sortedSizes{2} < edgeValues(iHistBin+1));

            % Downsample WT if there are more WT values than KO
            if length(WTindex) > length(KOindex)
                randomWTindex = randperm(length(WTindex), length(KOindex));
                DSsize_WT = [DSsize_WT, sortedSizes{1}(WTindex(randomWTindex))];
                DSslope_WT = [DSslope_WT, sortedSlopes{1}(WTindex(randomWTindex))];
                DSsize_KO = [DSsize_KO, sortedSizes{2}(KOindex)];
                DSslope_KO = [DSslope_KO, sortedSlopes{2}(KOindex)];
            % Downsample KO if there are more KO values than WT
            elseif length(KOindex) > length(WTindex)
                randomKOindex = randperm(length(KOindex), length(WTindex));
                DSsize_KO = [DSsize_KO, sortedSizes{2}(KOindex(randomKOindex))];
                DSslope_KO = [DSslope_KO, sortedSlopes{2}(KOindex(randomKOindex))];
                DSsize_WT = [DSsize_WT, sortedSizes{1}(WTindex)];
                DSslope_WT = [DSslope_WT, sortedSlopes{1}(WTindex)];
            % If the sizes are equal, just append all values
            else
                DSsize_WT = [DSsize_WT, sortedSizes{1}(WTindex)];
                DSslope_WT = [DSslope_WT, sortedSlopes{1}(WTindex)];
                DSsize_KO = [DSsize_KO, sortedSizes{2}(KOindex)];
                DSslope_KO = [DSslope_KO, sortedSlopes{2}(KOindex)];
            end
        end

        % Perform KS test for slopes and sizes
        [~, p_slope(iRandomIteration)] = kstest2(DSslope_WT, DSslope_KO);
        [~, p_size(iRandomIteration)] = kstest2(DSsize_WT, DSsize_KO);

    end
    data.populationData(1).phasePrecession.control.matchDistribution.byMean.dsSizes = DSsize_WT;
    data.populationData(2).phasePrecession.control.matchDistribution.byMean.dsSizes = DSsize_KO;
    data.populationData(1).phasePrecession.control.matchDistribution.byMean.dsSlopes = DSslope_WT;
    data.populationData(2).phasePrecession.control.matchDistribution.byMean.dsSlopes = DSslope_KO;
    data.populationData(1).phasePrecession.control.matchDistribution.byMean.sizePvalues = p_size;
    data.populationData(1).phasePrecession.control.matchDistribution.byMean.slopePvalues = p_slope;
    
    %% Step 4: Linear Regression Model 1 - Genotype Only 
    slopes_wt = data.populationData(1).phasePrecession.Slopes;
    slopes_ko = data.populationData(2).phasePrecession.Slopes;
    phase_slopes = [slopes_wt, slopes_ko]; 
    wt_genotype = ones(1,length(slopes_wt)); 
    ko_genotype = 2*ones(1,length(slopes_ko)); 
    Genotype = [wt_genotype, ko_genotype]; 
    validIdx = ~isnan(phase_slopes) & ~isnan(Genotype);
    phase_slopes = phase_slopes(validIdx); 
    Genotype = Genotype(validIdx); 

    tbl = table(phase_slopes', Genotype', ...
        'VariableNames', {'Slope', 'Genotype'});
    mdl = fitlm(tbl, 'Slope ~ Genotype');
    data.populationData(1).phasePrecession.control.linearRegression.genotypeOnly = mdl;
    
    %% Step 5: Linear Regression Model 2 - Field Size Only
    fieldSizes_wt = data.populationData(1).phasePrecession.MedianFieldSizes;
    fieldSizes_ko = data.populationData(2).phasePrecession.MedianFieldSizes;
    slopes_wt = data.populationData(1).phasePrecession.Slopes;
    slopes_ko = data.populationData(2).phasePrecession.Slopes;
    validIdx = ~isnan(fieldSizes_wt) & ~isnan(slopes_wt);
    field_sizes_wt = fieldSizes_wt(validIdx);
    phase_slopes_wt = slopes_wt(validIdx);
    validIdx = ~isnan(fieldSizes_ko) & ~isnan(slopes_ko);
    field_sizes_ko = fieldSizes_ko(validIdx);
    phase_slopes_ko = slopes_ko(validIdx);

    % Log-transformed regression
    log_FieldSize = [log(field_sizes_wt), log(field_sizes_ko)];
    phase_slopes = [phase_slopes_wt, phase_slopes_ko]; 
    tbl = table(phase_slopes', log_FieldSize', 'VariableNames', {'Slope', 'LogFieldSize'});
    mdl = fitlm(tbl, 'Slope ~ LogFieldSize');
    data.populationData(1).phasePrecession.control.linearRegression.fieldSizeOnly = mdl;
    
    %% Step 6: Linear Regression Model 3 - Theta MVL only
    theta_wt = data.populationData(1).MVL; 
    theta_ko = data.populationData(2).MVL;
    slopes_wt = data.populationData(1).phasePrecession.Slopes;
    slopes_ko = data.populationData(2).phasePrecession.Slopes;
    validIdx = ~isnan(theta_wt) & ~isnan(slopes_wt);
    theta_wt = theta_wt(validIdx);
    phase_slopes_wt = slopes_wt(validIdx);
    validIdx = ~isnan(theta_ko) & ~isnan(slopes_ko);
    theta_ko = theta_ko(validIdx);
    phase_slopes_ko = slopes_ko(validIdx);

    mvl = [theta_wt, theta_ko]; 
    phase_slopes = [phase_slopes_wt, phase_slopes_ko]; 
    tbl = table(phase_slopes', mvl', 'VariableNames', {'Slope', 'MVL'});
    mdl = fitlm(tbl, 'Slope ~ MVL');
    data.populationData(1).phasePrecession.control.linearRegression.MVLonly = mdl;
    
    %% Step 7: Linear Regression Model 4 - Full Model
    slopes_wt = data.populationData(1).phasePrecession.Slopes;
    slopes_ko = data.populationData(2).phasePrecession.Slopes;
    phase_slopes = [slopes_wt, slopes_ko]; 
    wt_genotype = ones(1,length(slopes_wt)); 
    ko_genotype = 2*ones(1,length(slopes_ko)); 
    Genotype = [wt_genotype, ko_genotype]; 
    fieldSizes_wt = data.populationData(1).phasePrecession.MedianFieldSizes;
    fieldSizes_ko = data.populationData(2).phasePrecession.MedianFieldSizes;
    fieldSizes = [fieldSizes_wt,fieldSizes_ko]; 
    theta_wt = data.populationData(1).MVL; 
    theta_ko = data.populationData(2).MVL;
    theta = [theta_wt, theta_ko]; 

    validIdx = ~isnan(phase_slopes) & ~isnan(Genotype) & ~isnan(fieldSizes) & ~isnan(theta);
    phase_slopes = phase_slopes(validIdx); 
    Genotype = Genotype(validIdx); 
    log_FieldSize = log(fieldSizes(validIdx)); 
    theta = theta(validIdx); 

    tbl = table(phase_slopes', Genotype', log_FieldSize', theta', ...
        'VariableNames', {'Slope', 'Genotype', 'FieldSize', 'Theta'});
    mdl = fitlm(tbl, 'Slope ~ Genotype + FieldSize + Theta');
    data.populationData(1).phasePrecession.control.linearRegression.fullModel.mdl = mdl;

    % VIFs for full model
    mdl_fs = fitlm(tbl, 'FieldSize ~ Genotype + Theta');
    R2_fs = mdl_fs.Rsquared.Ordinary;
    vif_FieldSize = 1 / (1 - R2_fs);
    data.populationData(1).phasePrecession.control.linearRegression.fullModel.vif_FieldSize = vif_FieldSize;

    mdl_theta = fitlm(tbl, 'Theta ~ Genotype + FieldSize');
    R2_theta = mdl_theta.Rsquared.Ordinary;
    vif_Theta = 1 / (1 - R2_theta);
    data.populationData(1).phasePrecession.control.linearRegression.fullModel.vif_theta = vif_Theta;

    mdl_geno = fitlm(tbl, 'Genotype ~ FieldSize + Theta');
    R2_geno = mdl_geno.Rsquared.Ordinary;
    vif_Genotype = 1 / (1 - R2_geno);
    data.populationData(1).phasePrecession.control.linearRegression.fullModel.vif_genotype = vif_Genotype;
    
    % Looking at interaction terms
    tbl = table(phase_slopes', Genotype', log_FieldSize', theta', ...
    'VariableNames', {'Slope', 'Genotype', 'FieldSize', 'ThetaMod'});
    mdl = fitlm(tbl, 'Slope ~ Genotype * FieldSize + Genotype * ThetaMod');
    data.populationData(1).phasePrecession.control.linearRegression.fullModel.mdl_interactions = mdl;
    
    %% Step 8: Looking at relationship between slope and field size using linear regression
    
    fieldSizes_wt = data.populationData(1).phasePrecession.MedianFieldSizes;
    fieldSizes_ko = data.populationData(2).phasePrecession.MedianFieldSizes;
    slopes_wt = data.populationData(1).phasePrecession.Slopes;
    slopes_ko = data.populationData(2).phasePrecession.Slopes;
    validIdx = ~isnan(fieldSizes_wt) & ~isnan(slopes_wt);
    field_sizes_wt = fieldSizes_wt(validIdx);
    phase_slopes_wt = slopes_wt(validIdx);
    validIdx = ~isnan(fieldSizes_ko) & ~isnan(slopes_ko);
    field_sizes_ko = fieldSizes_ko(validIdx);
    phase_slopes_ko = slopes_ko(validIdx);

    % Log-transformed regression
    wt_genotype = ones(1,length(field_sizes_wt)); 
    ko_genotype = 2*ones(1,length(field_sizes_ko)); 
    Genotype = [wt_genotype, ko_genotype]; 
    log_FieldSize = [log(field_sizes_wt), log(field_sizes_ko)];
    phase_slopes = [phase_slopes_wt, phase_slopes_ko]; 
    tbl = table(phase_slopes', log_FieldSize', Genotype', 'VariableNames', {'Slope', 'LogFieldSize', 'Genotype'});
    mdl = fitlm(tbl, 'Slope ~ LogFieldSize + Genotype');
    data.populationData(1).phasePrecession.control.linearRegression = mdl;
    
    %% Step 9: Bootstrap the median slopes to determine whether there's enough cells
    numIterations = settings.phasePrecession.control.bootstrapN;
    for iIteration = 1:numIterations;
        for iGenotype = 1:length(data.populationData)
            % Bootstrap the data
            withoutNaN = data.populationData(iGenotype).phasePrecession.Slopes;
            withoutNaN = withoutNaN(~isnan(withoutNaN));
            resampedData{iGenotype} = randsample(withoutNaN, length(withoutNaN), true); 
        end
        [~, bootstrappedPValue(iIteration)] = kstest2(resampedData{1}, resampedData{2}); 
        bootstrappedDiff(iIteration) = nanmean(resampedData{1}) - nanmean(resampedData{2}); 
    end
    
    data.populationData(1).phasePrecession.control.bootstrapping.Diff = bootstrappedDiff;
    data.populationData(1).phasePrecession.control.bootstrapping.pValue = bootstrappedPValue;
    
    %% Step 10: Shuffle the phases and the positions and get the slopes
    for iGenotype = 1:length(fieldnames(data.cellData));
        genotypes = fieldnames(data.cellData); 
        genotypeData = data.cellData.(genotypes{iGenotype}); 
        % Run analysis for high-firing cells only
        FRdata = genotypeData.highFiring;
        shPos_allIterations = []; shPhs_allIterations = [];
        numIterations = 1;
        for iIteration = 1:numIterations; 
            display(['On iteration ', num2str(iIteration), ' out of ', num2str(numIterations)]); 
            shPhs_populationSlopes = []; shPos_populationSlopes = []; 
            for iAnimal = 1:length(FRdata); 
                % Skip if empty
                if isempty(FRdata{iAnimal}) == 1; 
                    continue
                else
                    [~,n] = size(FRdata{iAnimal});
                    for iCluster = 1:n;
                        % Skip if empty
                        if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                            continue
                        else

                            % Assign variables based on running direction
                            directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original);
                            for iDir = 1:length(directions);
                                % Extract variables based on running direction
                                outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir), 'phasePrecession');
                                spkPhs = outputData.spkPhs; spkPos = outputData.spkPos; binnedSpkPos = outputData.binnedSpkPos;
                                spikeTimes = outputData.spikesByDirection;

                                % Loop through fields
                                numFieldsToAnalyze = whichField(settings.phasePrecession.fieldsToAnalyze, spkPhs);
                                trialPvalue = {}; 
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
                                                    max(abs(diff(spikeTimes{iField}{iTrial}))) <= settings.phasePrecession.ISIthreshold && ...
                                                    max(spikeTimes{iField}{iTrial}) - min(spikeTimes{iField}{iTrial}) >= settings.phasePrecession.timeRange;

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

                                                % Now shuffle the spkPhs and the spkPos and calculate the phase precession
                                                inputData.spkPhs = spkPhs{iField}{iTrial}; 
                                                inputData.spkPos = spkPosInput{iField}{iTrial}(randperm(length(spkPosInput{iField}{iTrial}))); 
                                                shPos_outputData = calculatePhasePrecession_v1_20241020(inputData, settings); 

                                                inputData.shSpkPhs = spkPhs{iField}{iTrial}(randperm(length(spkPhs{iField}{iTrial}))); 
                                                inputData.spkPos = spkPosInput{iField}{iTrial};
                                                shPhs_outputData = calculatePhasePrecession_v1_20241020(inputData, settings); 

                                                shPos_trialSlope{iField}(iTrial) = shPos_outputData.trialSlope; 
                                                shPhs_trialSlope{iField}(iTrial) = shPhs_outputData.trialSlope; 
                                                shPos_trialPvalue{iField}(iTrial) = shPos_outputData.p; 
                                                shPhs_trialPvalue{iField}(iTrial) = shPhs_outputData.p; 

                                                % If the slope is below the significance threshold, save the data
                                                if shPos_trialPvalue{iField}(iTrial) > settings.phasePrecession.significanceThreshold; 
                                                    shPos_trialSlope{iField}(iTrial) = NaN;
                                                end
                                                if shPhs_trialPvalue{iField}(iTrial) > settings.phasePrecession.significanceThreshold; 
                                                    shPhs_trialSlope{iField}(iTrial) = NaN;
                                                end
                                            else
                                                shPos_trialSlope{iField}(iTrial) = NaN;
                                                shPhs_trialSlope{iField}(iTrial) = NaN;
                                            end
                                        end
                                    end

                                    % If there are enough trials with slopes calculated, get the median of all slopes
                                    if sum(~isnan(shPos_trialSlope{iField})) >= settings.phasePrecession.trialThreshold; 
                                        shPos_slopeMedian(iField) = nanmedian(shPos_trialSlope{iField});
                                        shPhs_slopeMedian(iField) = nanmedian(shPhs_trialSlope{iField});
                                    else 
                                        shPos_slopeMedian(iField) = NaN; 
                                        shPhs_slopeMedian(iField) = NaN; 
                                    end
                                end

                                % Organize the population data and analysis
                                shPos_populationSlopes = [shPos_populationSlopes, shPos_slopeMedian];
                                shPhs_populationSlopes = [shPhs_populationSlopes, shPhs_slopeMedian];
                            end
                        end
                    end
                end
            end
            shPos_allIterations = [shPos_allIterations; shPos_populationSlopes]; 
            shPhs_allIterations = [shPhs_allIterations; shPhs_populationSlopes]; 
        end
        data.populationData(iGenotype).phasePrecession.control.shuffles.shuffPos = shPos_allIterations;
        data.populationData(iGenotype).phasePrecession.control.shuffles.shuffPhs = shPhs_allIterations;
    end
    
    
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
