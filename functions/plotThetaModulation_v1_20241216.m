function data = plotThetaModulation_v1_20241216(data, settings)
    % Generates plots associated with theta modulation
    % Written by Anja Payne
    % Last Modified: 02/28/2025

    % Inputs:
    %   1) data: the matlab structure where the theta modulation data is
    %      stored
    %   2) settings: the settings file where all theta-related settings are
    %      stored
    
    % Outputs:
    %   figure 1: the LFP, theta-filtered LFP, and spikes
    %   figure 2: the rose plots showing phase preference of all spikes
    %       plotted for each individual neuron
    %   figure 3: scatter plot on a rose plot showing the preferred phase
    %       and theta modulation of each neuron, for WT
    %   figure 4: scatter plot on a rose plot showing the preferred phase
    %       and theta modulation of each neuron, for KO
    %   figure 5: the preferred phase for KO and WT populations, displayed
    %       as a scatter plot around two concentric unit circles, with
    %       stats
    %   figure 6: the CDF of the mean vector length with stats for WT and
    %       KO
    %   figure 7: the MVL across the normalized field for WT and KO
    %       populations with a subplot of the p-values across the spatial
    %       bins
    %   figure 8: the CDF of the in-field spike numbers with stats for WT
    %       and KO populations
    %   figure 9: the CDF of the MVL for the bursts and singles, with
    %       stats, for the WT and KO populations
    
    % Steps/Figures:
    %   0) Get the folder the user wants to save plots into and ask them to
    %      select the figures they want to generate
    %   1) Plot the LFP with spikes for either all data or the animal/cell
    %      the user specified in the settings
    %   2) Plot the phase of the spikes for each cell as a rose plot
    %      histogram 
    %   3) Plot the scatter plot rose plots for WT and KO populations
    %   4) Plot the preferred phase for WT and KO populations
    %   5) Plot the CDF of the mean vector length for WT and KO populations
    %   6) Plot the mean +/- SEM MVL across the normalized field for WT and
    %      KO populations
    %   7) Plot the CDF of the in-field spike numbers for WT and KO
    %   8) Plot the CDF of the MVL for the bursts and singles for WT and KO
    %      populations
    
   
    close all; 
    
    if strcmp(settings.theta.plot.display, 'yes') == 1; 
        % Get the folders to save plots into
        mainFolder = uigetdir('C:\', 'Please select the folder you would like theta plots saved into.');
        figureSettings = getFigureFolders(mainFolder); 
        
        % Have the user select which plots they want to generate
        listOfFigures = getFiguresToPlot();
        
        %% Step 1: LFP with spikes
        if ismember(1, listOfFigures) == 1; 
            % Get list of genotypes to plot over
            if strcmp(settings.theta.plot.genotypes, 'all') == 1; 
                genotypeList = 1:length(fieldnames(data.cellData));
            else
                genotypeList = settings.theta.plot.genotypes; 
            end
            for iGenotype = genotypeList;
                genotypes = fieldnames(data.cellData); 
                genotypeData = data.cellData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;
                
                % Get animals to plot over
                if strcmp(settings.theta.plot.animals, 'all') == 1; 
                    animalList = 1:length(FRdata); 
                else
                    animalList = settings.theta.plot.animals; 
                end
                
                for iAnimal = animalList; 
                    % Skip if empty
                    if isempty(FRdata{iAnimal}) == 1; 
                        continue
                    else
                        
                        % Get cells to plot over
                        if strcmp(settings.theta.plot.cells, 'all') == 1; 
                            [~,n] = size(FRdata{iAnimal});
                            cellList = 1:n; 
                        else
                            cellList = settings.theta.plot.cells; 
                        end
                        
                        for iCluster = cellList; 
                            % Skip if empty
                            if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                                display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                                continue
                            else
                                
                                % Get directions to plot over
                                directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original);
                                if strcmp(settings.theta.plot.direction, 'all') == 1; 
                                    directionList = 1:length(directions); 
                                else
                                    directionList = settings.theta.plot.direction; 
                                end
                                
                                for iDir = directionList;
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                                    
                                    % Loop through fields
                                    if strcmp(settings.theta.fieldsToAnalyze, 'all fields') == 1;
                                        numField = length(FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.cw); 
                                    elseif strcmp(settings.theta.fieldsToAnalyze, 'best field') == 1;
                                        numField = 1; 
                                    end
                                    
                                    for iField = 1:numField;
                                    
                                        % Load LFP data
                                        load(FRdata{iAnimal}(iCluster).theta.LFP); 

                                        % Get the variables for plotting
                                        lfp_times = thetaLFP.LFPtimes;
                                        lfp = thetaLFP.LFPsamples;
                                        signal_filtered = thetaLFP.filteredLFP;
                                        if strcmp(directions(iDir), 'cw') == 1; 
                                            allTrials = FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.cw{iField}; 
                                            tSp = []; 
                                            for i = 1:length(allTrials); 
                                                tSp = [tSp; allTrials{i}]; 
                                            end
                                        elseif strcmp(directions(iDir), 'ccw') == 1; 
                                            allTrials = FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.ccw{1}; 
                                            tSp = []; 
                                            for i = 1:length(allTrials); 
                                                tSp = [tSp; allTrials{i}]; 
                                            end
                                        end

                                        % Plot!
                                        figures.LFP = figure(1); clf; hold on;
                                        plot(lfp_times/1000, lfp, 'k')
                                        plot(lfp_times/1000, signal_filtered, 'b'); 
                                        scatter(tSp/1000, 0.5*ones(1, length(tSp)), 200, '.r');
                                        xlabel('Time (sec)'); 
                                        ylabel('Arbitrary Units'); 
                                        set(gca, 'FontSize', 12); 
                                        set(figures.LFP, 'Position', [100, 200, 1800, 800]);

                                        % Save the figure
                                        figureSettings.filePath = figureSettings.filePath.LFP;
                                        figureSettings.name = [genotypes{iGenotype}, '_Animal', num2str(iAnimal), '_Cluster', ...
                                            num2str(iCluster), '_', directions{iDir}, '_Field', num2str(iField)];
                                        figureSettings.appendedFolder.binary = 'yes'; 
                                        figureSettings.appendedFolder.name = figureSettings.fileNameBase.LFP;
                                        figureSettings.fileTypes = {'fig', 'tiff'};
                                        saveFigure_v1_20240902(figures.LFP, figureSettings)
                                
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        %% Step 2: Individual rose plots

        if ismember(2, listOfFigures) == 1; 
            % Get list of genotypes to plot over
            if strcmp(settings.theta.plot.genotypes, 'all') == 1; 
                genotypeList = 1:length(fieldnames(data.cellData));
            else
                genotypeList = settings.theta.plot.genotypes; 
            end
            for iGenotype = genotypeList;
                genotypes = fieldnames(data.cellData); 
                genotypeData = data.cellData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;

                % Get animals to plot over
                if strcmp(settings.theta.plot.animals, 'all') == 1; 
                    animalList = 1:length(FRdata); 
                else
                    animalList = settings.theta.plot.animals; 
                end

                for iAnimal = animalList; 
                    % Skip if empty
                    if isempty(FRdata{iAnimal}) == 1; 
                        continue
                    else

                        % Get cells to plot over
                        if strcmp(settings.theta.plot.cells, 'all') == 1; 
                            [~,n] = size(FRdata{iAnimal});
                            cellList = 1:n; 
                        else
                            cellList = settings.theta.plot.cells; 
                        end

                        for iCluster = cellList; 
                            % Skip if empty
                            if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                                display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                                continue
                            else

                                % Get directions to plot over
                                directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original);
                                if strcmp(settings.theta.plot.direction, 'all') == 1; 
                                    directionList = 1:length(directions); 
                                else
                                    directionList = settings.theta.plot.direction; 
                                end

                                for iDir = directionList;
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);

                                    % Loop through fields
                                    if strcmp(settings.theta.fieldsToAnalyze, 'all fields') == 1;
                                        numField = length(FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.cw); 
                                    elseif strcmp(settings.theta.fieldsToAnalyze, 'best field') == 1;
                                        numField = 1; 
                                    end

                                    for iField = 1:numField;
                                        % Get variables to plot
                                        allPhases = []; 
                                        for i = 1:length(FRdata{iAnimal}(iCluster).theta.allSpikes.phases.cw{iField})
                                            allPhases = [allPhases, FRdata{iAnimal}(iCluster).theta.allSpikes.phases.cw{iField}{i}];
                                        end
                                        
                                        % Plot!
                                        figures.individualRose = figure(2); clf; 
                                        rose((allPhases), 24)
                                        set(gca, 'FontSize', 12); 
                                        set(figures.individualRose, 'Position', [100, 200, 500, 500]);
                                       
                                        % Save the figure
                                        figureSettings.filePath = figureSettings.filePath.individualRose;
                                        figureSettings.name = [genotypes{iGenotype}, '_Animal', num2str(iAnimal), '_Cluster', ...
                                            num2str(iCluster), '_', directions{iDir}, '_Field', num2str(iField)];
                                        figureSettings.appendedFolder.binary = 'yes'; 
                                        figureSettings.appendedFolder.name = figureSettings.fileNameBase.individualRose;
                                        figureSettings.fileTypes = {'fig', 'tiff'};
                                        saveFigure_v1_20240902(figures.individualRose, figureSettings)

                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        %% Step 3: Rose plot scatter plot for all neurons
        
        if ismember(3, listOfFigures) == 1; 
            figureSettings.filePath = figureSettings.filePath.population;
            for iGenotype = 1:length(fieldnames(data.cellData));
                genotypes = fieldnames(data.cellData); 

                % Plot scatter plot for all cells in population
                preferred_phase = deg2rad(data.populationData(iGenotype).preferredDirection); 
                mean_vector_length = data.populationData(iGenotype).MVL;

                if iGenotype == 1; 
                    figures.populationRose = figure(3); clf;
                    pointAppearance = 'ko';
                elseif iGenotype == 2; 
                    figures.populationRose = figure(4); clf; 
                    pointAppearance = 'go'; 
                end;
                polar(0, 0.8, 'w.'); % Forces the axis by placing a white dot
                hold on;        
                for i = 1:length(preferred_phase)
                    % Use polar function to plot each point
                    polar(preferred_phase(i), mean_vector_length(i), pointAppearance);
                end
                
                % Save the figure
                figureSettings.name = [genotypes{iGenotype}];
                figureSettings.appendedFolder.binary = 'no'; 
                figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
                figureSettings.fileTypes = {'fig', 'tiff'};
                saveFigure_v1_20240902(figures.populationRose, figureSettings);
            end
        end
        
        %% Step 4: Preferred Phase
        
        if ismember(4, listOfFigures) == 1; 
            figureSettings.filePath = figureSettings.filePath.population;
            
            % Assign the preferredPhase to WT and KO
            preferredPhase_WT = deg2rad(data.populationData(1).preferredDirection); 
            preferredPhase_KO = deg2rad(data.populationData(2).preferredDirection); 

            % Convert angles to Cartesian coordinates
            x_WT = cos(preferredPhase_WT);
            y_WT = sin(preferredPhase_WT);
            x_KO = 0.75.*cos(preferredPhase_KO); % Scaled down to nest inside WT
            y_KO = 0.75.*sin(preferredPhase_KO); % Scaled down to nest inside WT

            % Scatter plot
            figures.populationPrefPhase = figure(5); clf; hold on;
            scatter(x_WT, y_WT, 50, 'Marker', 'o', 'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerFaceColor', [0.5, 0.5, 0.5]); 
            scatter(x_KO, y_KO, 50, 'Marker', 'o', 'MarkerEdgeColor', [0.0, 0.5, 0.0], 'MarkerFaceColor', [0.0, 0.5, 0.0]); 
            axis equal; % Equal scaling for both axes

            % Add the mean +/- SEM as larger, filled circles
            WTmean = rad2deg(circ_mean(deg2rad(preferredPhase_WT'))); % Use circular stats toolbox since this is circular data
            KOmean = rad2deg(circ_mean(deg2rad(preferredPhase_KO'))); % Use circular stats toolbox since this is circular data
            WT_SEM = rad2deg(circ_std(deg2rad(preferredPhase_WT')) / sqrt(length(preferredPhase_WT)));
            KO_SEM = rad2deg(circ_std(deg2rad(preferredPhase_KO')) / sqrt(length(preferredPhase_KO)));

            WTstats = [WTmean-WT_SEM, WTmean, WTmean+WT_SEM];
            KOstats = [KOmean-KO_SEM, KOmean, KOmean+KO_SEM];
            xStats_WT = cos(WTstats); yStats_WT = sin(WTstats);
            xStats_KO = 0.75.*cos(KOstats); yStats_KO = 0.75.*sin(KOstats);
            
            scatter(xStats_WT, yStats_WT, 1500, 'Marker', '.', 'MarkerEdgeColor', [0.2, 0.2, 0.2], 'MarkerFaceColor', [0.2, 0.2, 0.2]); 
            scatter(xStats_KO, yStats_KO, 1500, 'Marker', '.', 'MarkerEdgeColor', [0.0, 0.2, 0.0], 'MarkerFaceColor', [0.0, 0.2, 0.0]); 

            % Display statistical significance
            % First, test for uniformity
            [pUnif_WT,~] = circ_rtest(preferredPhase_WT')
            [pUnif_KO,~] = circ_rtest(preferredPhase_KO') 
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not uniform, possibility it fits von Mises but check visually'); 
            end
            % Second, check that the concentration parameters are the same
            k_WT = 1 / circ_std(deg2rad(preferredPhase_WT')); 
            k_KO = 1 / circ_std(deg2rad(preferredPhase_KO'));  
            display('Check that the concentration parameters are similar'); 
            display(['Concentration parameter k for WT: ', num2str(k_WT)]);
            display(['Concentration parameter k for KO: ', num2str(k_KO)]);
            % If data passes assumptions, test for significance
            [p,~] = circ_wwtest(deg2rad(preferredPhase_WT'), deg2rad(preferredPhase_KO'));
            pValueDisplay = ['P-value is ', num2str(p), ' using circ wwtest'];
            display(pValueDisplay);
            display('Note that circ_wwtest assumes von Mises distribution and similar concentration parameters'); 
            
            % Save the figure
            figureSettings.name = 'PreferredPhase_WTandKO';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);
            saveFigure_v1_20240902(figures.populationPrefPhase, figureSettings);
        end
        
        %% Step 5: CDF of mean vector length
        
        if ismember(5, listOfFigures) == 1; 
            figureSettings.filePath = figureSettings.filePath.population;
            
            % Assign the preferredPhase to WT and KO
            MVL_WT = data.populationData(1).MVL; 
            MVL_KO = data.populationData(2).MVL; 

            % Plot
            figures.populationMVLcdf = figure(6); clf; hold on;
            plt_WT = cdfplot(MVL_WT); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(MVL_KO); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('MVL'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(MVL_WT); 
            KOmean = nanmean(MVL_KO);
            WT_SEM = nanstd(MVL_WT)/sqrt(length(MVL_WT)); 
            KO_SEM = nanstd(MVL_KO)/sqrt(length(MVL_KO)); 
            WT_y = interp1(sort(MVL_WT), linspace(0, 1, length(MVL_WT)), WTmean, 'linear', 'extrap');
            KO_y = interp1(sort(MVL_KO), linspace(0, 1, length(MVL_KO)), KOmean, 'linear', 'extrap');
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(MVL_WT)
            [h, pUnif_KO] = adtest(MVL_KO)
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(MVL_WT, MVL_KO);
            pValueDisplay = ['P-value is ', num2str(p), ' using kstest'];
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.name = 'MVL_CDF_WTandKO';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationMVLcdf, figureSettings);
        end
   
        %% Step 6: MVL across normalized field
        if ismember(6, listOfFigures) == 1; 
            figureSettings.filePath = figureSettings.filePath.population;
            
            for iGenotype = 1:length(fieldnames(data.cellData));
                normalizedMVL = []; normalizedMVL_firstHalf = []; normalizedMVL_secondHalf = []; 

                genotypes = fieldnames(data.cellData); 
                genotypeData = data.cellData.(genotypes{iGenotype}); 
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
                                
                                directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original);
                                for iDir = 1:length(directions);
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);

                                    % Get variables to plot
                                    allPhases = FRdata{iAnimal}(iCluster).theta.phases.(directions{iDir}); 
                                    bins = FRdata{iAnimal}(iCluster).inField.inFieldBinnedSpkPos.(directions{iDir}); 
                                    rateMap = FRdata{iAnimal}(iCluster).rateMaps.trialAverageMap.(directions{iDir}); 
                                    barcode = FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original.(directions{iDir}); 

                                    % Loop through fields
                                    if strcmp(settings.theta.fieldsToAnalyze, 'all fields') == 1;
                                        numField = length(FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.cw); 
                                    elseif strcmp(settings.theta.fieldsToAnalyze, 'best field') == 1;
                                        numField = 1; 
                                    end

                                    for iField = 1:numField;
                                        % Normalize the rate map (so it goes from 0 to 1 in
                                        % X and Y)
                                        rateMap = rateMap(barcode == iField); 
                                        rateMap = rateMap/nanmax(rateMap); 
                                        
                                        % Go through all the trials and sort the spikes
                                        % into bins
                                        byBinPhases = cell(1,67); 
                                        if length(allPhases{iField}) > 28;
                                            for iTrial = 1:length(allPhases{iField});
                                                for iSpikes = 1:length(allPhases{iField}{iTrial});
                                                    byBinPhases{bins{iField}{iTrial}(iSpikes)} = [byBinPhases{bins{iField}{iTrial}(iSpikes)}; allPhases{iField}{iTrial}(iSpikes)];
                                                end
                                            end
                                        end
                                        
                                        % Go through all the bins and get the theta
                                        % modulation
                                        mvlByBin = NaN(1,67); 
                                        for iBin = 1:length(byBinPhases)
                                            if isempty(byBinPhases{iBin}) == 0; 
                                                if length(byBinPhases{iBin}) ~= 1; 
                                                    mvlByBin(iBin) = circ_r(byBinPhases{iBin});
                                                end
                                            end
                                        end
                                        mvlByBinInField = mvlByBin(barcode == iField);
                                        [~, peakIndex] = nanmax(rateMap);
                                        firstHalf_mvl = mvlByBinInField(1:peakIndex);
                                        secondHalf_mvl = mvlByBinInField(peakIndex:end);
                                        
                                        % Interpolate all MVL values so they go from 0
                                        % to 67
                                        if length(mvlByBinInField) > 1 && length(firstHalf_mvl) > 1 && length(secondHalf_mvl) > 1; 
                                            newX_queryPoints = linspace(1, length(mvlByBinInField), 50);  
                                            newX = interp1([1:length(mvlByBinInField)], mvlByBinInField, newX_queryPoints);
                                            normalizedMVL = [normalizedMVL; newX]; 

                                            newX_queryPoints = linspace(1, length(firstHalf_mvl), 50); 
                                            newX_firstHalf = interp1([1:length(firstHalf_mvl)], firstHalf_mvl, newX_queryPoints);
                                            normalizedMVL_firstHalf = [normalizedMVL_firstHalf; newX_firstHalf];

                                            newX_queryPoints = linspace(1, length(secondHalf_mvl), 50); 
                                            newX_secondHalf = interp1([1:length(secondHalf_mvl)], secondHalf_mvl, newX_queryPoints);
                                            normalizedMVL_secondHalf = [normalizedMVL_secondHalf; newX_secondHalf]; 
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                normalizedMVL_byGenotype{iGenotype} = normalizedMVL;
                normalizedMVL_centeredAtPeak_byGenotype{iGenotype} = [normalizedMVL_firstHalf, normalizedMVL_secondHalf];
                
            end
            
            WTmean = nanmean(normalizedMVL_centeredAtPeak_byGenotype{1});
            WTsem = nanstd(normalizedMVL_centeredAtPeak_byGenotype{1})/sqrt(length(normalizedMVL_centeredAtPeak_byGenotype{1}));
            KOmean = nanmean(normalizedMVL_centeredAtPeak_byGenotype{2});
            KOsem = nanstd(normalizedMVL_centeredAtPeak_byGenotype{2})/sqrt(length(normalizedMVL_centeredAtPeak_byGenotype{2}));

            figures.MVLvsLocation = figure(7); clf; subplot(4,1,1:3); hold on; 
            plot([1/100:1/100:1], nanmean(normalizedMVL_centeredAtPeak_byGenotype{1}), 'k', 'LineWidth', 2); 
            plot([1/100:1/100:1], nanmean(normalizedMVL_centeredAtPeak_byGenotype{2}), 'Color', [0, 0.5, 0], 'LineWidth', 2);

            % Get SEM for WT
            x = [1/100:1/100:1];
            curve1 = WTmean + WTsem;
            curve2 = WTmean - WTsem;
            x2 = [x, fliplr(x)];
            inBetween = [curve1, fliplr(curve2)];
            fill(x2, inBetween, 'k', 'FaceAlpha', 0.3);

            % Get SEM for KO
            curve1 = KOmean + KOsem;
            curve2 = KOmean - KOsem;
            x2 = [x, fliplr(x)];
            inBetween = [curve1, fliplr(curve2)];
            fill(x2, inBetween, 'g', 'FaceAlpha', 0.3);
            xlim([0,1]);
            ylabel('Mean Vector Length');
            set(gca, 'FontSize', 14); 
            
            pValue = [];
            for i = 1:100; 
                WT = normalizedMVL_centeredAtPeak_byGenotype{1}(:,i); 
                KO = normalizedMVL_centeredAtPeak_byGenotype{2}(:,i); 
                [~, pValue(i)] = kstest2(WT, KO);
            end
            subplot(4,1,4);
            semilogy([1/100:1/100:1], pValue, 'r'); hold on; plot([1/100, 1], [0.05, 0.05], '--r');
            xlim([0,1]);
            ylabel('P-value'); xlabel('Normalized Field Location'); 
            set(gca, 'FontSize', 14); 
            
            % Save the figure
            figureSettings.name = 'MVLvsLocation_WTandKO';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.MVLvsLocation, figureSettings);
        end
        
        %% Step 7: CDF of spike number
        if ismember(7, listOfFigures) == 1; 
            figureSettings.filePath = figureSettings.filePath.population;
            
            for iGenotype = 1:length(fieldnames(data.cellData));
                genotypes = fieldnames(data.cellData); 
                genotypeData = data.cellData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;
                spikeNumbers = []; 

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
                                directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original);
                                for iDir = 1:length(directions);
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                                    allSpikes_allTrials = []; 
                                    
                                    % Loop through fields
                                    if strcmp(settings.theta.fieldsToAnalyze, 'all fields') == 1;
                                        numField = length(FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.cw); 
                                    elseif strcmp(settings.theta.fieldsToAnalyze, 'best field') == 1;
                                        numField = 1; 
                                    end

                                    for iField = 1:numField;
                                        % Get variables to plot
                                        spikesByDirection = FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.(directions{iDir}){iField}; 
                                        
                                        for iTrial = 1:length(spikesByDirection); 
                                           allSpikes_allTrials = [allSpikes_allTrials; spikesByDirection{iTrial}];  
                                        end
                                        spikeNumbers = [spikeNumbers, sum(~isnan(allSpikes_allTrials))];
                                    end
                                end
                            end
                        end
                    end
                end
                all_spikeNumbers{iGenotype} = spikeNumbers; 
            end
            
            % Plot and save for all spikes
            figures.populationSpkNumbersCDF = figure(8); clf; hold on;
            plt_WT = cdfplot(all_spikeNumbers{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(all_spikeNumbers{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Number of In-Field Spikes'); 
            set(gca, 'FontSize', 14); 
            WTmean = nanmean(all_spikeNumbers{1}); 
            KOmean = nanmean(all_spikeNumbers{2});
            WT_SEM = nanstd(all_spikeNumbers{1})/sqrt(length(all_spikeNumbers{1})); 
            KO_SEM = nanstd(all_spikeNumbers{2})/sqrt(length(all_spikeNumbers{2}));
            [f, x] = ecdf(all_spikeNumbers{1});
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(all_spikeNumbers{2});
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx); 
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);
            [h, pUnif_WT] = adtest(all_spikeNumbers{1})
            [h, pUnif_KO] = adtest(all_spikeNumbers{2})
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            [~, p] = kstest2(all_spikeNumbers{1}, all_spikeNumbers{2});
            pValueDisplay = ['P-value is ', num2str(p), ' using kstest'];
            display(pValueDisplay);
            annotation('textbox', [0.5, 0.15, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);
            % Save the figure
            figureSettings.name = 'NumberOfSpikes_CDF_WTandKO';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationSpkNumbersCDF, figureSettings);
        end
        
        %% Step 8: CDF of mean vector length for bursts and singles
        
        if ismember(8, listOfFigures) == 1; 
            figureSettings.filePath = figureSettings.filePath.population;
            
            % Assign the preferredPhase to WT and KO
            MVL_WT_bursts = data.populationData(1).burstsMVL; 
            MVL_KO_bursts = data.populationData(2).burstsMVL; 
            MVL_WT_singles = data.populationData(1).singlesMVL; 
            MVL_KO_singles = data.populationData(2).singlesMVL; 

            % Plot bursts
            figures.populationMVLcdfBursts = figure(9); clf; hold on;
            plt_WT = cdfplot(MVL_WT_bursts); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(MVL_KO_bursts); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('MVL'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(MVL_WT_bursts); 
            KOmean = nanmean(MVL_KO_bursts);
            WT_SEM = nanstd(MVL_WT_bursts)/sqrt(length(MVL_WT_bursts)); 
            KO_SEM = nanstd(MVL_KO_bursts)/sqrt(length(MVL_KO_bursts)); 
            WT_y = interp1(sort(MVL_WT_bursts), linspace(0, 1, length(MVL_WT_bursts)), WTmean, 'linear', 'extrap');
            KO_y = interp1(sort(MVL_KO_bursts), linspace(0, 1, length(MVL_KO_bursts)), KOmean, 'linear', 'extrap');
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(MVL_WT_bursts)
            [h, pUnif_KO] = adtest(MVL_KO_bursts)
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(MVL_WT_bursts, MVL_KO_bursts);
            pValueDisplay = ['P-value is ', num2str(p), ' using kstest'];
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.name = 'MVLbursts_CDF_WTandKO';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationMVLcdfBursts, figureSettings);
            
            % Plot singles
            figures.populationMVLcdfSingles = figure(10); clf; hold on;
            plt_WT = cdfplot(MVL_WT_singles); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(MVL_KO_singles); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('MVL'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(MVL_WT_singles); 
            KOmean = nanmean(MVL_KO_singles);
            WT_SEM = nanstd(MVL_WT_singles)/sqrt(length(MVL_WT_singles)); 
            KO_SEM = nanstd(MVL_KO_singles)/sqrt(length(MVL_KO_singles)); 
            WT_y = interp1(sort(MVL_WT_singles), linspace(0, 1, length(MVL_WT_singles)), WTmean, 'linear', 'extrap');
            KO_y = interp1(sort(MVL_KO_singles), linspace(0, 1, length(MVL_KO_singles)), KOmean, 'linear', 'extrap');
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(MVL_WT_singles)
            [h, pUnif_KO] = adtest(MVL_KO_singles)
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(MVL_WT_singles, MVL_KO_singles);
            pValueDisplay = ['P-value is ', num2str(p), ' using kstest'];
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.name = 'MVLsingles_CDF_WTandKO';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationMVLcdfSingles, figureSettings);
        end
        
    end
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function figureSettings = getFigureFolders(mainFolder)

    % For the LFP with spikes
    figureSettings.fileNameBase.LFP = 'thetaLFPwithSpikes';
    figureSettings.filePath.LFP = getMostRecentFilePath_v1_20240723(figureSettings.fileNameBase.LFP, '', mainFolder);
    
    % For the individual rose plots
    figureSettings.fileNameBase.individualRose = 'thetaRosePlots_individualCells';
    figureSettings.filePath.individualRose = getMostRecentFilePath_v1_20240723(figureSettings.fileNameBase.individualRose, '', mainFolder);

    % For the population rose plots
    figureSettings.fileNameBase.population = 'theta_populationPlots';
    figureSettings.filePath.population = getMostRecentFilePath_v1_20240723(figureSettings.fileNameBase.population, '', mainFolder);

end

function selectedPlots = getFiguresToPlot()
   % Define the list of available plots
    plotOptions = {'LFP with theta-filtered LFP and spikes', ...
        'rose plots of spiking for each cell', ...
        'rose plot scatter plot of preferred phase vs. MVL for all cells', ...
        'preferred phase for WT and KO populations'...
        'CDF of mean vector length for WT and KO populations', ...
        'MVL across normalized field'...
        'CDF of number of spikes for WT and KO populations', ...
        'CDF of MVL for bursts and singles'}; 
    
    % Display a dialog box to select the plots
    selectedPlots = listdlg('ListString', plotOptions, ...
    	'SelectionMode', 'multiple', ...
    	'PromptString', 'Select the plots you want to generate:', ...
    	'ListSize', [300, 100]);
end
    
    
    