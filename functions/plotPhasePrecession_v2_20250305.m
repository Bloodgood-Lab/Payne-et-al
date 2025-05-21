function plotPhasePrecession_v2_20250305(data, settings)
    % Generates plots related to phase precession analysis
    % Written by Anja Payne
    % Last Modified: 03/11/2025
    
    % Inputs:
    %   1) data: the matlab structure where the in-field theta phases
    %      by-trial are saved
    %   2) settings: the settings file where the settings for the phase
    %      precession calculations are saved

    % Outputs:
    %   figure 1: the LFP, theta-filtered LFP, spikes, and phase preference
    %   figure 2: 
    %     
    
    % Steps/Figures:
    %   0) Get the folder the user wants to save plots into and ask them to
    %      select the figures they want to generate
    %   1) Plot the LFP with spikes and theta phases for either all data or
    %      the animal/cell the user specified in the settings
    %   2) Plot the spikes with slopes for each trial as well as the
    %      histogram of the slopes
    %   3) Plot the CDF of the median slopes
    %   4) Get the scatter plots of the median field size and the median
    %      field slope
    %   5) Compare the size of the place fields associated with the trial 
    %      with the median slope 
    
    close all;
    
    if strcmp(settings.phasePrecession.plot.display, 'yes') == 1;
        % Get the folders to save plots into
        mainFolder = uigetdir('C:\', 'Please select the folder you would like phase precession plots saved into.');
        figureSettings = getFigureFolders(mainFolder);
    
        % Have the user select which plots they want to generate
        listOfFigures = getFiguresToPlot();
    
        %% Step 1: LFP with spikes and theta phase
        if ismember(1, listOfFigures) == 1; 
            % Get list of genotypes to plot over
            if strcmp(settings.phasePrecession.plot.genotypes, 'all') == 1; 
                genotypeList = 1:length(fieldnames(data.cellData));
            else
                genotypeList = settings.phasePrecession.plot.genotypes; 
            end
            for iGenotype = genotypeList;
                genotypes = fieldnames(data.cellData); 
                genotypeData = data.cellData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;
                
                % Get animals to plot over
                if strcmp(settings.phasePrecession.plot.animals, 'all') == 1; 
                    animalList = 1:length(FRdata); 
                else
                    animalList = settings.phasePrecession.plot.animals; 
                end
                
                for iAnimal = animalList; 
                    % Skip if empty
                    if isempty(FRdata{iAnimal}) == 1; 
                        continue
                    else
                        
                        % Get cells to plot over
                        if strcmp(settings.phasePrecession.plot.cells, 'all') == 1; 
                            [~,n] = size(FRdata{iAnimal});
                            cellList = 1:n; 
                        else
                            cellList = settings.phasePrecession.plot.cells; 
                        end
                        
                        for iCluster = cellList; 
                            % Skip if empty
                            if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                                display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                                continue
                            else
                                
                                % Get directions to plot over
                                directions = fieldnames(FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes);
                                if strcmp(settings.phasePrecession.plot.direction, 'all') == 1; 
                                    directionList = 1:length(directions); 
                                else
                                    directionList = settings.phasePrecession.plot.direction; 
                                end
                                
                                for iDir = directionList;
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                                    
                                    % Loop through fields
                                    if strcmp(settings.phasePrecession.fieldsToAnalyze, 'all fields') == 1;
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
                                            allTrials_spkTimes = FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.cw{iField}; 
                                            allTrials_spkPhases = FRdata{iAnimal}(iCluster).theta.allSpikes.phases.cw{iField}; 
                                            tSp = []; spkPhs = []; 
                                            for i = 1:length(allTrials_spkTimes); 
                                                tSp = [tSp; allTrials_spkTimes{i}]; 
                                                spkPhs = [spkPhs, allTrials_spkPhases{i}]; 
                                            end
                                        elseif strcmp(directions(iDir), 'ccw') == 1; 
                                            allTrials_spkTimes = FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.ccw{1}; 
                                            allTrials_spkPhases = FRdata{iAnimal}(iCluster).theta.phases.ccw{iField}; 
                                            tSp = []; spkPhs = []; 
                                            for i = 1:length(allTrials_spkTimes); 
                                                tSp = [tSp; allTrials_spkTimes{i}]; 
                                                spkPhs = [spkPhs, allTrials_spkPhases{i}]; 
                                            end
                                        end

                                        % Plot!
                                        figures.LFP = figure(1); clf; 
                                        subplot(2,1,1); hold on;
                                        plot(lfp_times/1000, lfp, 'k')
                                        plot(lfp_times/1000, signal_filtered, 'b'); 
                                        scatter(tSp/1000, 0.5*ones(1, length(tSp)), 200, '.r');
                                        xlabel('Time (sec)'); 
                                        ylabel('Arbitrary Units'); 
                                        set(gca, 'FontSize', 12); 
                                        set(figures.LFP, 'Position', [100, 200, 1800, 800]);
                                        subplot(2,1,2); 
                                        scatter(tSp, spkPhs, '.r');

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
        
        %% Step 2: Spikes with slopes for each trial; histogram of slopes 
        %  across trials
        if ismember(2, listOfFigures) == 1; 
            % Get list of genotypes to plot over
            if strcmp(settings.phasePrecession.plot.genotypes, 'all') == 1; 
                genotypeList = 1:length(fieldnames(data.cellData));
            else
                genotypeList = settings.phasePrecession.plot.genotypes; 
            end
            for iGenotype = genotypeList;
                genotypes = fieldnames(data.cellData); 
                genotypeData = data.cellData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;
                
                % Get animals to plot over
                if strcmp(settings.phasePrecession.plot.animals, 'all') == 1; 
                    animalList = 1:length(FRdata); 
                else
                    animalList = settings.phasePrecession.plot.animals; 
                end
                
                for iAnimal = animalList; 
                    % Skip if empty
                    if isempty(FRdata{iAnimal}) == 1; 
                        continue
                    else
                        
                        % Get cells to plot over
                        if strcmp(settings.phasePrecession.plot.cells, 'all') == 1; 
                            [~,n] = size(FRdata{iAnimal});
                            cellList = 1:n; 
                        else
                            cellList = settings.phasePrecession.plot.cells; 
                        end
                        
                        for iCluster = cellList; 
                            % Skip if empty
                            if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                                display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                                continue
                            else
                                
                                % Get directions to plot over
                                directions = fieldnames(FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes);
                                if strcmp(settings.phasePrecession.plot.direction, 'all') == 1; 
                                    directionList = 1:length(directions); 
                                else
                                    directionList = settings.phasePrecession.plot.direction; 
                                end
                                
                                for iDir = directionList;
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                                    
                                    % Extract the needed data from structure
                                    outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir), 'plotPhasePrecession');
                                    outputData.genotype = genotypes{iGenotype}; outputData.animal = iAnimal; outputData.cell = iCluster; outputData.dir = directions{iDir}; 
                                    outputData.figureSettings = figureSettings;

                                    % Plot 
                                    plotSpikesAndSlopes(outputData, settings); 
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
        
        %% Step 3: Plot the CDF of the median slopes
        if ismember(3, listOfFigures) == 1; 
            medianSlopes_WT = data.populationData(1).phasePrecessionSlopes; 
            medianSlopes_KO = data.populationData(2).phasePrecessionSlopes; 

            % Plot
            figures.populationMedianSlopeCDF = figure(3); clf; hold on;
            plt_WT = cdfplot(medianSlopes_WT); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(medianSlopes_KO); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('MVL'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(medianSlopes_WT); 
            KOmean = nanmean(medianSlopes_KO);
            WT_SEM = nanstd(medianSlopes_WT)/sqrt(length(medianSlopes_WT)); 
            KO_SEM = nanstd(medianSlopes_KO)/sqrt(length(medianSlopes_KO)); 
            [f, x] = ecdf(medianSlopes_WT);
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(medianSlopes_KO);
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(medianSlopes_WT)
            [h, pUnif_KO] = adtest(medianSlopes_KO)
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(medianSlopes_WT, medianSlopes_KO);
            pValueDisplay = ['P-value is ', num2str(p), ' using kstest'];
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'MedianSlopes_CDF_WTandKO';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationMedianSlopeCDF, figureSettings);
        end
        
        %% Step 4: Get scatter plots of the median field size and the median field slope
        if ismember(4, listOfFigures) == 1; 
            for iGenotype = 1:length(fieldnames(data.cellData));
                
                if iGenotype == 1; 
                    figures.populationSlopeVsSize = figure(4); hold on; %clf;
                    colorAppearance = [0.5 0.5 0.5];
                    figureSettings.name = 'sizeVsSlope_WT';
                    figureSettings.filePath = figureSettings.filePath.population;
                    figureSettings.appendedFolder.binary = 'no'; 
                    figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
                    figureSettings.fileTypes = {'fig', 'tiff'};
                elseif iGenotype == 2; 
                    figures.populationSlopeVsSize = figure(4); %clf; 
                    colorAppearance = [0.0 0.5 0.0]; 
                    figureSettings.name = 'sizeVsSlope_KO';
                end;
            
            medianSize = data.populationData(iGenotype).phasePrecessionMedianFieldSizes;
            medianSlope = data.populationData(iGenotype).phasePrecessionSlopes; 
            scatter(medianSize, rad2deg(medianSlope), 75, colorAppearance, 'filled'); hold on; 
            xlabel('Median Place Field Size (cm)'); 
            ylabel('Median Slope (degrees/cm)');
            set(gca, 'FontSize', 14); 

            % Save the figure
            saveFigure_v1_20240902(figures.populationSlopeVsSize, figureSettings);
            end

        end
        
        %% Step 5: Compare the median field size
        if ismember(5, listOfFigures) == 1; 
            medianSizes_WT = data.populationData(1).phasePrecessionMedianFieldSizes; 
            medianSizes_KO = data.populationData(2).phasePrecessionMedianFieldSizes; 

            % Plot
            figures.populationMedianFieldSizeeCDF = figure(6); clf; hold on;
            plt_WT = cdfplot(medianSizes_WT); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(medianSizes_KO); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('MVL'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(medianSizes_WT); 
            KOmean = nanmean(medianSizes_KO);
            WT_SEM = nanstd(medianSizes_WT)/sqrt(length(medianSizes_WT)); 
            KO_SEM = nanstd(medianSizes_KO)/sqrt(length(medianSizes_KO)); 
            [f, x] = ecdf(medianSizes_WT);
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(medianSizes_KO);
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(medianSizes_WT)
            [h, pUnif_KO] = adtest(medianSizes_KO)
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(medianSizes_WT, medianSizes_KO);
            pValueDisplay = ['P-value is ', num2str(p), ' using kstest'];
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'MedianFieldSize_CDF_WTandKO';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationMedianFieldSizeeCDF, figureSettings);
        end
        
        %% Step 6: Compare the median field size
        if ismember(6, listOfFigures) == 1; 
            medianSizes_WT = data.populationData(1).phasePrecessionMedianFieldSizes; 
            medianSizes_KO = data.populationData(2).phasePrecessionMedianFieldSizes; 

            % Plot
            figures.populationMedianFieldSizeeCDF = figure(7); clf; hold on;
            plt_WT = cdfplot(medianSizes_WT); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(medianSizes_KO); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Median Field Size (cm)'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(medianSizes_WT); 
            KOmean = nanmean(medianSizes_KO);
            WT_SEM = nanstd(medianSizes_WT)/sqrt(length(medianSizes_WT)); 
            KO_SEM = nanstd(medianSizes_KO)/sqrt(length(medianSizes_KO)); 
            [f, x] = ecdf(medianSizes_WT);
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(medianSizes_KO);
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(medianSizes_WT)
            [h, pUnif_KO] = adtest(medianSizes_KO)
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(medianSizes_WT, medianSizes_KO);
            pValueDisplay = ['P-value is ', num2str(p), ' using kstest'];
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'MedianFieldSize_CDF_WTandKO';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationMedianFieldSizeeCDF, figureSettings);
        end
        
        %% Step 7: Compare the average field size
        if ismember(7, listOfFigures) == 1; 
            averageSizes_WT = data.populationData(1).phasePrecession.averageFieldSizes; 
            averageSizes_KO = data.populationData(2).phasePrecession.averageFieldSizes; 
            slopes_WT = data.populationData(1).phasePrecessionSlopes; 
            slopes_KO = data.populationData(2).phasePrecessionSlopes; 
            if length(averageSizes_WT) ~= length(slopes_WT)
                display('Error, unequal number of fields in WT slope and field arrays'); 
            elseif length(averageSizes_KO) ~= length(slopes_KO)
                display('Error, unequal number of fields in KO slope and field arrays'); 
            else
                averageSizes_WT(isnan(slopes_WT)) = [];
                averageSizes_KO(isnan(slopes_KO)) = [];
            end
            length(averageSizes_WT)
            length(averageSizes_KO)
            
            % Plot
            figures.populationAverageFieldSizeeCDF = figure(8); clf; hold on;
            plt_WT = cdfplot(averageSizes_WT); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(averageSizes_KO); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Trial-Averaged Field Size (cm)'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(averageSizes_WT); 
            KOmean = nanmean(averageSizes_KO);
            WT_SEM = nanstd(averageSizes_WT)/sqrt(length(averageSizes_WT)); 
            KO_SEM = nanstd(averageSizes_KO)/sqrt(length(averageSizes_KO)); 
            [f, x] = ecdf(averageSizes_WT);
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(averageSizes_KO);
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(averageSizes_WT)
            [h, pUnif_KO] = adtest(averageSizes_KO)
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(averageSizes_WT, averageSizes_KO);
            pValueDisplay = ['p=', num2str(p), ' using kstest'];
            display(pValueDisplay);
            annotation('textbox', [0.60, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'TrialAveragedFieldSize_CDF_WTandKO';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationAverageFieldSizeeCDF, figureSettings);
        end
        
        %% Step 8: Match the distributions of field size and compare slopes
        if ismember(8, listOfFigures) == 1; 
            WTsizeDistribution = data.populationData(1).phasePrecession.control.matchDistribution.byMedian.dsSizes; 
            KOsizeDistribution = data.populationData(2).phasePrecession.control.matchDistribution.byMedian.dsSizes; 
            WTslopeDistribution = data.populationData(1).phasePrecession.control.matchDistribution.byMedian.dsSlopes; 
            KOslopeDistribution = data.populationData(2).phasePrecession.control.matchDistribution.byMedian.dsSlopes; 
            slopePvalues = data.populationData(1).phasePrecession.control.matchDistribution.byMedian.slopePvalues;
            
            % Plot the distributions of the sizes
            figures.controlMedianSizeDistribution = figure(9); clf; hold on;
            histogram(WTsizeDistribution, [0:4:100], 'FaceColor', [0.5 0.5 0.5]); 
            histogram(KOsizeDistribution, [0:4:100], 'FaceColor', [0.0 0.5 0.0]); 
            xlabel('Downsampled Median Field Sizes'); 
            ylabel('Counts')
            set(gca, 'FontSize', 14); 

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'Control_SizeDistribution';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.controlMedianSizeDistribution, figureSettings);
            
            % Plot the distributions of the sizes
            figures.controlMedianSizeSlopeDistribution = figure(10); clf; hold on;
            plt_WT = cdfplot(WTslopeDistribution); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(KOslopeDistribution); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Downsampled Slopes (degrees/cm)'); 
            set(gca, 'FontSize', 14); 
            
            % Save the figure
            figureSettings.name = 'Control_SlopeDistribution';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.controlMedianSizeSlopeDistribution, figureSettings);
            
            % Plot the CDF of the pValues
            figures.controlMedianSizeDistribution_pValues = figure(11); clf; hold on;
            plt = cdfplot(slopePvalues); 
            set(plt, 'Color', [0 0 0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('P-Values between Slopes'); 
            set(gca, 'FontSize', 14); 
            
            % Save the figure
            figureSettings.name = 'Control_pValues';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.controlMedianSizeDistribution_pValues, figureSettings);
        end
        
        %% Step 9: Match the distributions of field size and compare slopes
        if ismember(9, listOfFigures) == 1; 
            %WTsizeDistribution = data.populationData(1).phasePrecession.control.matchDistribution.byMean.dsSizes; 
            %KOsizeDistribution = data.populationData(2).phasePrecession.control.matchDistribution.byMean.dsSizes; 
            WTsizeDistribution = data.populationData(1).phasePrecessionMedianFieldSizes; 
            KOsizeDistribution = data.populationData(2).phasePrecessionMedianFieldSizes;
            WTslopeDistribution = data.populationData(1).phasePrecession.control.matchDistribution.byMean.dsSlopes; 
            KOslopeDistribution = data.populationData(2).phasePrecession.control.matchDistribution.byMean.dsSlopes; 
            slopePvalues = data.populationData(1).phasePrecession.control.matchDistribution.byMean.slopePvalues;
            
            % Plot the distributions of the sizes
            figures.controlMeanSizeDistribution = figure(12); clf; hold on;
            histogram(WTsizeDistribution, [0:4:100], 'FaceColor', [0.5 0.5 0.5]); 
            histogram(KOsizeDistribution, [0:4:100], 'FaceColor', [0.0 0.5 0.0]); 
            xlabel('Downsampled Mean Field Sizes'); 
            ylabel('Counts')
            set(gca, 'FontSize', 14); 

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'Control_SizeDistribution';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.controlMeanSizeDistribution, figureSettings);
            
            % Plot the distributions of the sizes
            figures.controlMeanSizeSlopeDistribution = figure(13); clf; hold on;
            plt_WT = cdfplot(WTslopeDistribution); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(KOslopeDistribution); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Downsampled Slopes (degrees/cm)'); 
            set(gca, 'FontSize', 14); 
            
            % Save the figure
            figureSettings.name = 'Control_SlopeDistribution';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.controlMeanSizeSlopeDistribution, figureSettings);
            
            % Plot the CDF of the pValues
            figures.controlMeanSizeDistribution_pValues = figure(14); clf; hold on;
            plt = cdfplot(slopePvalues); 
            set(plt, 'Color', [0 0 0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('P-Values between Slopes'); 
            set(gca, 'FontSize', 14); 
            
            % Save the figure
            figureSettings.name = 'Control_pValues';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.controlMeanSizeDistribution_pValues, figureSettings);
        end
        
        %% Step 10: Visualize how many cells are excluded with each criteria
        if ismember(10, listOfFigures) == 1; 
            WT_percentSpatBin = sum(data.populationData(1).phasePrecession.thresholds.spatialBin)/length(data.populationData(1).phasePrecession.thresholds.spatialBin); 
            KO_percentSpatBin = sum(data.populationData(2).phasePrecession.thresholds.spatialBin)/length(data.populationData(2).phasePrecession.thresholds.spatialBin); 
            WT_percentISI = sum(data.populationData(1).phasePrecession.thresholds.ISI)/length(data.populationData(1).phasePrecession.thresholds.ISI); 
            KO_percentISI = sum(data.populationData(2).phasePrecession.thresholds.ISI)/length(data.populationData(2).phasePrecession.thresholds.ISI); 
            WT_percenttimeRange = sum(data.populationData(1).phasePrecession.thresholds.timeRange)/length(data.populationData(1).phasePrecession.thresholds.timeRange); 
            KO_percenttimeRange = sum(data.populationData(2).phasePrecession.thresholds.timeRange)/length(data.populationData(2).phasePrecession.thresholds.timeRange); 
            WT_percentP = sum(data.populationData(1).phasePrecession.thresholds.pValue)/length(data.populationData(1).phasePrecession.thresholds.pValue); 
            KO_percentP = sum(data.populationData(2).phasePrecession.thresholds.pValue)/length(data.populationData(2).phasePrecession.thresholds.pValue); 
            WT_trialNum = sum(data.populationData(1).phasePrecession.thresholds.trialNum)/length(data.populationData(1).phasePrecession.thresholds.trialNum); 
            KO_trialNum = sum(data.populationData(2).phasePrecession.thresholds.trialNum)/length(data.populationData(2).phasePrecession.thresholds.trialNum); 
            allData = 100.*[WT_percentSpatBin, KO_percentSpatBin, ...
                WT_percentISI, KO_percentISI, ...
                WT_percenttimeRange, KO_percenttimeRange, ...
                WT_percentP, KO_percentP, ...
                WT_trialNum, KO_trialNum]; 
            
            % Plot the distributions of the sizes
            figures.exclusionCriteria = figure(15); clf; hold on;
            x = 1:length(allData); 
            for i = 1:length(allData); 
                if mod(i,2) == 1;
                    bar(x(i), allData(i), 'FaceColor', [0.5, 0.5, 0.5]); 
                else
                    bar(x(i), allData(i), 'FaceColor', [0.0, 0.5, 0.0]); 
                end
            end
            %xlabel('Downsampled Mean Field Sizes'); 
            ylabel('Percent of Cells Included')
            set(gca, 'XTick', [])  % Hide default tick labels
            text(1.5, -0.05*max(allData), 'spatialBin', 'HorizontalAlignment', 'center')
            text(3.5, -0.05*max(allData), 'ISI', 'HorizontalAlignment', 'center')
            text(5.5, -0.05*max(allData), 'timeRange', 'HorizontalAlignment', 'center')
            text(7.5, -0.05*max(allData), 'p-Value', 'HorizontalAlignment', 'center')
            text(9.5, -0.05*max(allData), 'trialNum', 'HorizontalAlignment', 'center')
            set(gca, 'FontSize', 14); 

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'Control_ExclusionCriteria';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.exclusionCriteria, figureSettings);
                  
        end
        
        %% Step 11: Report how many trials are included for each neuron
        if ismember(11, listOfFigures) == 1; 
            for iGenotype = 1:length(fieldnames(data.cellData));
                genotypes = fieldnames(data.cellData); 
                genotypeData = data.cellData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;
                numTrials{iGenotype} = []; 
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
                                directions = fieldnames(FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes);
                                for iDir = 1:length(directions);
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                                    
                                    % Loop through fields
                                    if strcmp(settings.phasePrecession.fieldsToAnalyze, 'all fields') == 1;
                                        numField = length(FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.(directions{iDir})); 
                                    elseif strcmp(settings.theta.fieldsToAnalyze, 'best field') == 1;
                                        numField = 1; 
                                    end
                                    
                                    for iField = 1:numField;
                                        % If the median slope is NaN, skip
                                        if isnan(FRdata{iAnimal}(iCluster).phasePrecession.medianSlope.(directions{iDir})); 
                                            continue;
                                        else
                                            cellNumTrials = sum(~isnan(FRdata{iAnimal}(iCluster).phasePrecession.allSlopes.(directions{iDir}){iField}));  
                                            numTrials{iGenotype} = [numTrials{iGenotype}, cellNumTrials]; 
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            % Plot the distributions of the number of trials
            figures.controlNumberTrials = figure(16); clf; hold on;
            plt_WT = cdfplot(numTrials{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(numTrials{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Number of Trials Included to get Mean'); 
            set(gca, 'FontSize', 14); 

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'Control_NumberOfTrials';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.controlNumberTrials, figureSettings);        
        end
        
        %% Step 12: Show standard deviations of slopes across trials
        if ismember(12, listOfFigures) == 1; 
            for iGenotype = 1:length(fieldnames(data.cellData));
                genotypes = fieldnames(data.cellData); 
                genotypeData = data.cellData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;
                standardDeviations{iGenotype} = [];
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
                                directions = fieldnames(FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes);
                                for iDir = 1:length(directions);
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                                    
                                    % Loop through fields
                                    if strcmp(settings.phasePrecession.fieldsToAnalyze, 'all fields') == 1;
                                        numField = length(FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.(directions{iDir})); 
                                    elseif strcmp(settings.theta.fieldsToAnalyze, 'best field') == 1;
                                        numField = 1; 
                                    end
                                    
                                    for iField = 1:numField;
                                        % If the median slope is NaN, skip
                                        if isnan(FRdata{iAnimal}(iCluster).phasePrecession.medianSlope.(directions{iDir})); 
                                            continue;
                                        else
                                            cellStandardDeviation = nanstd(FRdata{iAnimal}(iCluster).phasePrecession.allSlopes.(directions{iDir}){iField});
                                            standardDeviations{iGenotype} = [standardDeviations{iGenotype}, cellStandardDeviation]; 
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            % Plot the distributions of the number of trials
            figures.controlSTDev = figure(17); clf; hold on;
            plt_WT = cdfplot(standardDeviations{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(standardDeviations{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Standard Deviation of Trials For Each Cell'); 
            set(gca, 'FontSize', 14); 
            
            % If data is not normal, perform kstest
            [~, p] = kstest2(standardDeviations{1}, standardDeviations{2});
            pValueDisplay = ['P-value is ', num2str(p), ' using kstest'];
            display(['p-value is: ', num2str(pValueDisplay)]);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);
            display(['WT N = ', num2str(length(standardDeviations{1})), '; KO N = ', num2str(length(standardDeviations{2}))]);

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'Control_StandardDeviation';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.controlSTDev, figureSettings);        
        end
        
        %% Step 13: Show results of bootstrapping the data
        if ismember(13, listOfFigures) == 1; 
            
            % Plot the distributions of the difference between WT mean and
            % KO mean
            figures.bootstrappedDiff = figure(18); clf; hold on;
            plt = histogram(data.populationData(1).phasePrecession.control.bootstrapping.Diff);
            set(plt, 'FaceColor', [0.5 0.5 0.5]); 
            xlabel('WT Mean Slope - KO Mean Slope'); 
            set(gca, 'FontSize', 14); 

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'Control_BootstrappedDiff';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.bootstrappedDiff, figureSettings); 
            
            % Plot the distributions of the p-value between WT and KO
            figures.bootstrappedPValue = figure(19); clf; hold on;
            plt = histogram(data.populationData(1).phasePrecession.control.bootstrapping.pValue);
            set(plt, 'FaceColor', [0.5 0.5 0.5]); 
            xlabel('p-value of bootstrapped data'); 
            set(gca, 'FontSize', 14); 
            
            % Save the figure
            figureSettings.name = 'Control_BootstrappedPValue';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.bootstrappedPValue, figureSettings); 
        end
        
        %% Step 14: Show the results of the linear regression
        if ismember(14, listOfFigures) == 1; 
            mdl = data.populationData(1).phasePrecession.control.linearRegression;
            figures.linearRegression = figure(20);
            t = uitable(figures.linearRegression, 'Data', mdl.Coefficients{:,:}, ...
            'ColumnName', mdl.Coefficients.Properties.VariableNames, ...
            'RowName', mdl.Coefficients.Properties.RowNames, ...
            'Units', 'Normalized', ...
            'Position', [0 0 1 1]);

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'Control_LinearRegression.tif';
            [figureSettings.filePath{1}, figureSettings.name]
            saveas(figures.linearRegression, [figureSettings.filePath{1}, figureSettings.name])
        end
        
        %% Step 15: Get the CDFs for the shuffled datasets
        if ismember(15, listOfFigures) == 1; 
            
            %%%%% Shuffled phases: Difference in means + p-values %%%%%
            
            % Get the difference in means and the p-values
            [m,~] = size(data.populationData(1).phasePrecession.control.shuffles.shuffPhs);  
            for iIteration = 1:m; 
                WTdata = data.populationData(1).phasePrecession.control.shuffles.shuffPhs(iIteration,:); 
                KOdata = data.populationData(2).phasePrecession.control.shuffles.shuffPhs(iIteration,:); 
                average(iIteration) = nanmean(WTdata)-nanmean(KOdata); 
                [~, pValue(iIteration)] = kstest2(WTdata, KOdata); 
            end
            
            % Plot the distributions of the difference between WT mean and
            % KO mean
            figures.shuffPhsDiff = figure(21); clf; hold on;
            edges = [-0.06:0.01:0.04]; 
            plt = histogram(average, edges); 
            set(plt, 'FaceColor', [0.5 0.5 0.5]); 
            xlabel('WT Mean Slope - KO Mean Slope'); 
            set(gca, 'FontSize', 14); 

            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'Control_ShuffledPhs_Diff';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.shuffPhsDiff, figureSettings); 
            
            % Plot the distributions of the p-value between WT and KO
            figures.shuffPhsPValue = figure(22); clf; hold on;
            edges = [0:0.025:0.8]; 
            plt = histogram(pValue, edges);
            set(plt, 'FaceColor', [0.5 0.5 0.5]); 
            xlabel('p-value of shuffled phase data'); 
            set(gca, 'FontSize', 14); 
            
            % Save the figure
            figureSettings.name = 'Control_ShuffledPhs_pValues';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.shuffPhsPValue, figureSettings); 
            
            %%%%% Shuffled positions: Difference in means + p-values %%%%%
            
            % Get the difference in means and the p-values
            [m,~] = size(data.populationData(1).phasePrecession.control.shuffles.shuffPos);  
            for iIteration = 1:m; 
                WTdata = data.populationData(1).phasePrecession.control.shuffles.shuffPos(iIteration,:); 
                KOdata = data.populationData(2).phasePrecession.control.shuffles.shuffPos(iIteration,:); 
                average(iIteration) = nanmean(WTdata)-nanmean(KOdata); 
                [~, pValue(iIteration)] = kstest2(WTdata, KOdata); 
            end
            
            % Plot the distributions of the difference between WT mean and
            % KO mean
            figures.shuffPosDiff = figure(23); clf; hold on;
            edges = [-0.06:0.01:0.04]; 
            plt = histogram(average, edges); 
            set(plt, 'FaceColor', [0.5 0.5 0.5]); 
            xlabel('WT Mean Slope - KO Mean Slope'); 
            set(gca, 'FontSize', 14); 

            % Save the figure
            figureSettings.name = 'Control_ShuffledPos_Diff';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.shuffPosDiff, figureSettings); 
            
            % Plot the distributions of the p-value between WT and KO
            figures.shuffPosPValue = figure(24); clf; hold on;
            edges = [0:0.025:0.8]; 
            plt = histogram(pValue, edges);
            set(plt, 'FaceColor', [0.5 0.5 0.5]); 
            xlabel('p-value of shuffled positions'); 
            set(gca, 'FontSize', 14); 
            
            % Save the figure
            figureSettings.name = 'Control_ShuffledPos_pValues';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.shuffPosPValue, figureSettings); 
        end
        
        %% Step 16: Plot the log of field size against slope 
        if ismember(16, listOfFigures) == 1; 
            % Format data
            fieldSizes_wt = data.populationData(1).phasePrecession.MedianFieldSizes;
            fieldSizes_ko = data.populationData(2).phasePrecession.MedianFieldSizes;
            slopes_wt = data.populationData(1).phasePrecession.Slopes;
            slopes_ko = data.populationData(2).phasePrecession.Slopes;
            validIdx = ~isnan(fieldSizes_wt) & ~isnan(slopes_wt);
            field_sizes_wt = log(fieldSizes_wt(validIdx));
            slopes_wt = slopes_wt(validIdx);
            validIdx = ~isnan(fieldSizes_ko) & ~isnan(slopes_ko);
            field_sizes_ko = log(fieldSizes_ko(validIdx));
            slopes_ko = slopes_ko(validIdx);
            
            % Get the correlation information
            [rho, pval] = corrcoef([field_sizes_wt, field_sizes_ko], [slopes_wt, slopes_ko]);
            corrDisplay = ['Rho = ', num2str(rho(2)), ' and p-value = ', num2str(pval(2))]
            
            % Plot the distributions of the p-value between WT and KO
            figures.scatterSizeVsSlope = figure(25); clf; hold on;
            scatter(field_sizes_wt, slopes_wt, 50, [0.5, 0.5, 0.5], 'filled');
            scatter(field_sizes_ko, slopes_ko, 50, [0.0, 0.5, 0.0], 'filled');
            xlabel('Slope (radians/cm)'); 
            ylabel('Log(Field Size)'); 
            set(gca, 'FontSize', 14); 
            annotation('textbox', [0.02, 0, 0.2, 0.03], 'String', ...
                ['Data File: ', settings.dataSavePath], 'FitBoxToText', 'on', 'BackgroundColor', 'none', ...
                'EdgeColor', 'none');
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', corrDisplay, ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);
            
            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'Control_ScatterLogOfSizeVsSlope';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.scatterSizeVsSlope, figureSettings); 
        end
        
        %% Step 17: Plot the theta modulation MVL against slope 
        if ismember(17, listOfFigures) == 1; 
            % Format data
            theta_mvl_wt = data.populationData(1).MVL;
            theta_mvl_ko = data.populationData(2).MVL;
            slopes_wt = data.populationData(1).phasePrecession.Slopes;
            slopes_ko = data.populationData(2).phasePrecession.Slopes;
            validIdx = ~isnan(theta_mvl_wt) & ~isnan(slopes_wt);
            mvl_wt = theta_mvl_wt(validIdx);
            slopes_wt = slopes_wt(validIdx);
            validIdx = ~isnan(theta_mvl_ko) & ~isnan(slopes_ko);
            mvl_ko = theta_mvl_ko(validIdx);
            slopes_ko = slopes_ko(validIdx);
            
            % Get the correlation information
            [rho, pval] = corrcoef([mvl_wt, mvl_ko], [slopes_wt, slopes_ko]); 
            corrDisplay = ['Rho = ', num2str(rho(2)), ' and p-value = ', num2str(pval(2))];
            
            % Plot the distributions of the p-value between WT and KO
            figures.scatterMVLvsSlope = figure(25); clf; hold on;
            scatter(mvl_wt, slopes_wt, 50, [0.5, 0.5, 0.5], 'filled');
            scatter(mvl_ko, slopes_ko, 50, [0.0, 0.5, 0.0], 'filled');
            xlabel('Slope (radians/cm)'); 
            ylabel('MVL'); 
            set(gca, 'FontSize', 14); 
            annotation('textbox', [0.02, 0, 0.2, 0.03], 'String', ...
                ['Data File: ', settings.dataSavePath], 'FitBoxToText', 'on', 'BackgroundColor', 'none', ...
                'EdgeColor', 'none');
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', corrDisplay, ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);
            
            % Save the figure
            figureSettings.filePath = figureSettings.filePath.population;
            figureSettings.name = 'Control_ScatterMVLvsSlope';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureSettings.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.scatterMVLvsSlope, figureSettings); 
        end
        
    end
end   
         
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selectedPlots = getFiguresToPlot()
   % Define the list of available plots
    plotOptions = {'LFP with spikes and theta phase',...
        'spikes with slopes and histogram of slopes',...
        'CDF of median slopes for WT and KO populations',...
        'scatter plots of median slope vs median field size',...
        'CDF of median field size compared between WT and KO populations',...
        'CDF of the average field size compared between WT and KO populations',...
        'CDF of the trial-averaged field size compared between WT and KO populations', ...
        'control analysis: match the median size distributions between WT and KO', ...
        'control analysis: match the mean size distributions between WT and KO', ...
        'control analysis: visualize how many cells are excluded by criteria', ...
        'control analysis: report number of trials for each cell', ...
        'control analysis: report standard deviation across trials', ...
        'control analysis: bootstrap the data to determine if N is sufficient', ...
        'control analysis: linear regression results', ...
        'control analysis: shuffled position and phase', ...
        'control analysis: scatter plot of field size against slope', ... 
        'control analysis: scatter plot of theta MVL against slope'}; 
    
    % Display a dialog box to select the plots
    selectedPlots = listdlg('ListString', plotOptions, ...
    	'SelectionMode', 'multiple', ...
    	'PromptString', 'Select the plots you want to generate:', ...
    	'ListSize', [300, 200]);
end

function figureSettings = getFigureFolders(mainFolder)

    % For the LFP with spikes and theta phases
    figureSettings.fileNameBase.LFP = 'phasePrecession_thetaLFPwithSpikesAndSpikePhases';
    figureSettings.filePath.LFP = getMostRecentFilePath_v1_20240723(figureSettings.fileNameBase.LFP, '', mainFolder);
    
    % For the spikes and slopes plots
    figureSettings.fileNameBase.spikesAndSlopes = 'phasePrecession_spikesAndSlopes';
    figureSettings.filePath.spikesAndSlopes = getMostRecentFilePath_v1_20240723(figureSettings.fileNameBase.spikesAndSlopes, '', mainFolder);

    % For the population plots
    figureSettings.fileNameBase.population = 'phasePrecession_populationPlots';
    figureSettings.filePath.population = getMostRecentFilePath_v1_20240723(figureSettings.fileNameBase.population, '', mainFolder);

end

function plotSpikesAndSlopes(inputData, settings)
    % Ask the user to select the folder to save figures into
    spkPhs = inputData.spkPhsForPlot; spkPos = inputData.spkPosForPlot; binnedSpkPos = inputData.binnedSpkPos; 
    slopeMedian = inputData.slopeMedian; allSlopes = inputData.allSlopes; 
    offsets = inputData.offsets; rSquared = inputData.rSquared;

    % Loop through fields
    if strcmp(settings.phasePrecession.fieldsToAnalyze, 'all fields') == 1;
        numFieldsToAnalyze = length(spkPhs);
    elseif strcmp(settings.phasePrecession.fieldsToAnalyze, 'best field') == 1;
        numFieldsToAnalyze = 1;
    end
    for iField = 1:numFieldsToAnalyze;
        close all; 
        % Set the values for the subplots
        [subplotNumbRows, allLinePlotRange, histoPlotRange] = determineSubplotLocation(length(spkPhs{iField}));

        % Loop through all the trials
        allPlotSlopes = []; allYoffsets = []; maxX = []; 
        for iTrial = 1:length(spkPhs{iField});
            figures.allTrialsFig = figure(2); set(figures.allTrialsFig, 'Position', [100, 200, 1800, 800]);
            if iTrial == 1; clf(figures.allTrialsFig); else hold on; end;

            % Determine subplot to plot onto
            if iTrial <= 10; iScatterSub = iTrial; 
            elseif iTrial > 10 && iTrial <= 20; iScatterSub = iTrial+3; 
            elseif iTrial > 20 && iTrial <= 30; iScatterSub = iTrial+6;
            elseif iTrial > 30; iScatterSub = iTrial+9; 
            end

            % Plot scatter subplot
            subplot(subplotNumbRows, 13, iScatterSub); hold on;
            scatter(spkPos{iField}{iTrial}, spkPhs{iField}{iTrial}, 200, '.k');
            set(gca, 'FontSize', 12);
            title(num2str(iTrial));
            ylim([0, 2*pi]);

            % If there are enough spatial bins
            if nanmax(binnedSpkPos{iField}{iTrial})-nanmin(binnedSpkPos{iField}{iTrial}) < settings.phasePrecession.spatialBinThreshold;                                             continue;
                continue;
            else
                % If there were spikes in-field that trial
                if length(spkPhs{iField}{iTrial}) >= settings.phasePrecession.spikeThreshold; 
                    % Get the fit values
                    plotSlope = allSlopes{iField}(iTrial);
                    yOffset = offsets{iField}(iTrial);
                    allPlotSlopes = [allPlotSlopes, plotSlope]; 
                    allYoffsets = [allYoffsets, yOffset]; 
                    maxX = nanmax([maxX, nanmax(spkPos{iField}{iTrial})]);

                    % Plot that trial's scatter plot
                    xPlot = [0:max(spkPos{iField}{iTrial})-min(spkPos{iField}{iTrial})]
                    if strcmp(settings.phasePrecession.circularity, 'shift') == 1 ...
                            || strcmp(settings.phasePrecession.circularity, 'none') == 1;
                        plotSlope
                        yOffset
                        yPlot = (xPlot*plotSlope) + yOffset
                        plot(xPlot, yPlot, 'r'); 
                        display('lines should be there')
                    elseif strcmp(settings.phasePrecession.circularity, 'doubling') == 1;
                        if strcmp(settings.phasePrecession.fit, 'linear') == 1;
                            yPlot1 = (xPlot*(plotSlope) + yOffset) - pi;
                            yPlot2 = (xPlot*(plotSlope) + yOffset) + pi;
                        elseif strcmp(settings.phasePrecession.fit, 'circular') == 1;
                            yPlot1 = (xPlot*(plotSlope) + yOffset) + pi;
                            yPlot2 = (xPlot*(plotSlope) + yOffset) + 3*pi;
                        end
                        plot(xPlot, yPlot1, 'r');
                        plot(xPlot, yPlot2, 'r'); 
                    end
                    if length(xPlot)>1; xlim([min(xPlot), max(xPlot)]); end
                    ylim([0, 2*pi]); 

                    % Plot the overlay of all slopes
                    subplot(subplotNumbRows, 13, allLinePlotRange); hold on;
                    if strcmp(settings.phasePrecession.circularity, 'shift') == 1;
                        plot(xPlot, yPlot, 'k');
                    elseif strcmp(settings.phasePrecession.circularity, 'doubling') == 1;
                        plot(xPlot, yPlot1, 'k'); 
                    end
                end
            end
        end   
        % Clean up the scatter plots
        cleanUpPlot(settings, slopeMedian(iField), nanmean(allSlopes{iField}), nanmean(rSquared{iField}));

        % Clean up the slope overlay plot and
        % plot the average slope
        if ~isnan(slopeMedian(iField)) && ~isnan(nanmean(allSlopes{iField}));
            inputData.subplotNumbRows = subplotNumbRows;
            inputData.allLinePlotRange = allLinePlotRange;
            inputData.maxX = maxX; inputData.allPlotSlopes = allPlotSlopes; 
            inputData.yOffset = allYoffsets; 
            inputData.median = slopeMedian(iField);
            inputData.mean = nanmean(allSlopes{iField});
            cleanUpAndPlotOverlay(inputData);
        end

        % Plot a histogram of all slopes
        plotHistogram(subplotNumbRows, histoPlotRange, allSlopes{iField}, slopeMedian(iField));

        % Save the figure
        figureSettings.filePath = inputData.figureSettings.filePath.spikesAndSlopes;
        figureSettings.name = [inputData.genotype, '_Animal', num2str(inputData.animal), '_Cluster', ...
            num2str(inputData.cell), '_', inputData.dir, '_Field', num2str(iField)];
        figureSettings.appendedFolder.binary = 'yes'; 
        figureSettings.appendedFolder.name = inputData.figureSettings.fileNameBase.spikesAndSlopes;
        figureSettings.fileTypes = {'fig', 'tiff'};
        saveFigure_v1_20240902(figures.allTrialsFig, figureSettings)
    end
end

function [numRows, linePlotRange, histoRange] = determineSubplotLocation(dataLength)
    numRows = ceil(dataLength/10);
    if dataLength <= 10;
        linePlotRange = [12, 13];
    elseif dataLength > 10 && dataLength <= 20;
        linePlotRange = [12, 13];
        histoRange = [25, 26]; 
    elseif dataLength > 20 && dataLength <= 30;
        linePlotRange = [12, 26]; 
        histoRange = [38, 39]; 
    elseif dataLength > 30
        linePlotRange = [12, 26]; 
        histoRange = [38, 52]; 
    end
end

function cleanUpPlot(settings, slopeMedian, slopeMean, rSquared)
    han = axes('Position', [0.125 0.125 0.8 0.8], 'Visible', 'off');
    han.YLabel.Visible = 'on';
    ylabel(han, 'Theta Phase (degrees)', 'FontSize', 20);
    han = axes('Position', [0.065 0.10 0.9 0.9], 'Visible', 'off');
    han.XLabel.Visible = 'on';
    xlabel(han, 'Linear Position', 'FontSize', 20);
    annotation('textbox', [0.02, 0, 0.2, 0.03], 'String', ...
        ['Data File: ', settings.dataSavePath], 'FitBoxToText', 'on', 'BackgroundColor', 'none', ...
        'EdgeColor', 'none');
    annotation('textbox', [0.75, 0.98, 0.2, 0.03], 'String', ...
        ['Median Slope: ', num2str(slopeMedian)], 'FitBoxToText', 'on', ...
        'BackgroundColor', 'none', 'EdgeColor', 'none', 'FontSize', 12);
    annotation('textbox', [0.87, 0.98, 0.2, 0.03], 'String', ...
        ['Mean Slope: ', num2str(slopeMean)], 'FitBoxToText', 'on', ...
        'BackgroundColor', 'none', 'EdgeColor', 'none', 'FontSize', 12);
    annotation('textbox', [0.75, 0.95, 0.2, 0.03], 'String', ...
        ['R-squared: ', num2str(rSquared)], 'FitBoxToText', 'on', ...
        'BackgroundColor', 'none', 'EdgeColor', 'none', 'FontSize', 12);
    set(gcf, 'PaperUnits', 'Inches', 'PaperPositionMode', 'auto');
end

function cleanUpAndPlotOverlay(inputData)
    subplot(inputData.subplotNumbRows, 13, inputData.allLinePlotRange); hold on;
    xAverageLine = 0:inputData.maxX;
    inputData.allPlotSlopes(isinf(inputData.allPlotSlopes)==1) = NaN;
    inputData.yOffset(isinf(inputData.yOffset)==1) = NaN;
    yMedianLine = xAverageLine*inputData.median + nanmean(inputData.yOffset);
    yMeanLine = xAverageLine*inputData.mean + nanmean(inputData.yOffset);
    medianLine = plot(xAverageLine, yMedianLine, '-r', 'LineWidth', 2);
    meanLine = plot(xAverageLine, yMeanLine, '--b', 'LineWidth', 2);
    %legend([medianLine, meanLine], {'Median', 'Mean'}, 'Location', 'northwest');
    legend('boxoff');
    ylim([0, 2*pi]); 
    ylabel('Theta Phase (degrees)');
    xlabel('Linear Position');  
    set(gca,'FontSize', 12);
end

function plotHistogram(numRows, plotRange, slopes, medianSlope)
    subplot(numRows, 13, plotRange); hold on;
    slopeHist = histogram(slopes, 20, 'Normalization', 'probability');
    maxValue = max(slopeHist.Values);
    medianLine = plot([medianSlope, medianSlope], [0, maxValue], '-r', 'LineWidth', 2); 
    meanLine = plot([nanmean(slopes), nanmean(slopes)], [0, maxValue], '--b', 'LineWidth', 2); 
    legend([medianLine, meanLine], {'Median', 'Mean'}, 'Location', 'northwest');
    legend('boxoff');
    ylabel('Number of Trials');
    xlabel('Phase Precession Slopes');  
    set(gca,'FontSize', 12); 
end

function [barcode, PFSize, PFNumber] = getPlaceFields(map, thresholds)
    % Identifies the bins that are in-field
    % Inputs: 
    %   1) map: the collapsed map to be analyzed
    %   2) settings: the settings file which should contain lowThresh: the
    %      lower threshold for what is considered a place field and
    %      highThresh: the maximum firing rate that the field should
    %      contain
    % Output: 
    %   1) barcode: the barcode that tells you which bins belong to each
    %      field
    %   2) PFSize: the size of place fields
    %   3) PFNumber: the number of place fields
    % Steps: 
    %   1) Find the bins above the low threshold 
    %   2) Only include place fields if the max is above the high threshold
    %   3) Re-number the fields so that they are ordered by firing rate 
    %      (i.e. field 1 has the highest FR, then field 2, etc.)
    %   4) Get the number of place fields 

    %% Step 1: Find the bins above the low threshold 
    low_threshold = thresholds(1) * nanmax(map(:)); 
    barcode = bwlabel(map > low_threshold);

    %% Step 2: Only include fields where the max is above the high thresh
    high_threshold = thresholds(2) * nanmax(map(:));
    tempSize = []; FR = []; count = 1; 
    for i = 1:max(max(barcode)); 
        if max(map(barcode == i)) >= high_threshold;
            tempSize(count) = 4*length(find(barcode == i)); % Multiply by 4 since bin size is 4 cm.
            FR(count) = 16000 * nanmax(map(barcode == i)); % Get the max firing rate so that the fields can be ordered by FR
            barcode(barcode == i) = count; 
            count = count + 1;
        elseif max(map(barcode == i)) < high_threshold;  
            barcode(barcode == i) = 0; 
        end;
    end

    %% Step 3: Reorder the fields in order of firing rate
    [~, sortedInd] = sort(FR, 'descend'); % Get the order of the fields from highest FR to lowest
    barcode = 10*barcode; % Need this to be a larger value so that you can reorganize and reset values starting at 1
    PFSize = []; 
    for iField = 1:length(sortedInd); 
        barcode(barcode == 10*sortedInd(iField)) = iField;
        PFSize(iField) = tempSize(sortedInd(iField));
    end

    %% Step 4: Get the number of place fields
    PFNumber = length(PFSize);
end
                    