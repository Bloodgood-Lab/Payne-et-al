function data = controlPhasePrecession_v1_20250228(data, settings)
    % Generates plots related to phase precession analysis
    % Written by Anja Payne
    % Last Modified: 03/11/2025
    
    % Inputs:
    %   1) data: the matlab structure where the in-field theta phases
    %      by-trial are saved
    %   2) settings: the settings file where the settings for the phase
    %      precession calculations are saved

    % Outputs:
    %   
    %     
    
    % Steps/Figures:
    %   1) Get the average size of the field for each cell that was
    %      included in the phase precession slope analysis
    %   2) Match the distributions of the place field sizes then compare 
    %      the slopes
    
    
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
        data.populationData(iGenotype).phasePrecessionControl.averageFieldSizes = allFieldSizes;
    end
    
    %% Step 2: Match the distributions of the place field sizes then compare the slopes
    
    for iGenotype = 1:length(data.populationData)
        
        % Sort by size and apply the same sorting to the slope vector
        [sortedSizes{iGenotype}, sortedSize_idx] = sort(data.populationData(iGenotype).phasePrecessionMedianFieldSizes);
        sortedSlopes{iGenotype} = data.populationData(iGenotype).phasePrecessionSlopes(sortedSize_idx);
        
        % Bin the size vector up
        edgeValues = [0:4:max(sortedSizes{1})];
        [sizeHistogramValues(iGenotype, :), ~] = histcounts(sortedSizes{iGenotype}, edgeValues);
        
    end
    
    % Set random seed for reproducibility
    rng('shuffle');

    % Initialize variables
    numIterations = 1000;
    maxN = length(sizeHistogramValues);
    edgeValues = [0:4:max(sortedSizes{1})];

    % Preallocate results
    p_slope = zeros(numIterations, 1);
    p_size = zeros(numIterations, 1);

    % Loop through random iterations
    for iRandomIteration = 1:numIterations
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Debugging
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        figure(1); hold on; 
        plt_WT = cdfplot(DSslope_WT);
        set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1); 
        plt_KO = cdfplot(DSslope_KO); 
        set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1); 
        sum(~isnan(DSslope_WT))
        sum(~isnan(DSslope_KO))
        %}
    end
    data.populationData(1).phasePrecessionControl.matchDistribution.byMedian.dsSizes = DSsize_WT;
    data.populationData(2).phasePrecessionControl.matchDistribution.byMedian.dsSizes = DSsize_KO;
    data.populationData(1).phasePrecessionControl.matchDistribution.byMedian.dsSlopes = DSslope_WT;
    data.populationData(2).phasePrecessionControl.matchDistribution.byMedian.dsSlopes = DSslope_KO;
    data.populationData(1).phasePrecessionControl.matchDistribution.byMedian.sizePvalues = p_size;
    data.populationData(1).phasePrecessionControl.matchDistribution.byMedian.slopePvalues = p_slope;
    
    %% Step 3: Match the distributions of the average place field sizes then compare the slopes
    sizeHistogramValues = []; 
    for iGenotype = 1:length(data.populationData)
        
        % Sort by size and apply the same sorting to the slope vector
        [sortedSizes{iGenotype}, sortedSize_idx] = sort(data.populationData(iGenotype).phasePrecession.averageFieldSizes);
        sortedSlopes{iGenotype} = data.populationData(iGenotype).phasePrecessionSlopes(sortedSize_idx);
        
        % Bin the size vector up
        edgeValues = [0:4:max(sortedSizes{1})];
        [sizeHistogramValues(iGenotype, :), ~] = histcounts(sortedSizes{iGenotype}, edgeValues);
        
    end
    
    % Set random seed for reproducibility
    rng('shuffle');

    % Initialize variables
    numIterations = 1000;
    maxN = length(sizeHistogramValues);
    edgeValues = [0:4:max(sortedSizes{1})];

    % Preallocate results
    p_slope = zeros(numIterations, 1);
    p_size = zeros(numIterations, 1);

    % Loop through random iterations
    for iRandomIteration = 1:numIterations
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Debugging
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        figure(1); hold on; 
        plt_WT = cdfplot(DSslope_WT);
        set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1); 
        plt_KO = cdfplot(DSslope_KO); 
        set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1); 
        sum(~isnan(DSslope_WT))
        sum(~isnan(DSslope_KO))
        %}
        
    end
    data.populationData(1).phasePrecessionControl.matchDistribution.byMean.dsSizes = DSsize_WT;
    data.populationData(2).phasePrecessionControl.matchDistribution.byMean.dsSizes = DSsize_KO;
    data.populationData(1).phasePrecessionControl.matchDistribution.byMean.dsSlopes = DSslope_WT;
    data.populationData(2).phasePrecessionControl.matchDistribution.byMean.dsSlopes = DSslope_KO;
    data.populationData(1).phasePrecessionControl.matchDistribution.byMean.sizePvalues = p_size;
    data.populationData(1).phasePrecessionControl.matchDistribution.byMean.slopePvalues = p_slope;
    
    
end
    
    

