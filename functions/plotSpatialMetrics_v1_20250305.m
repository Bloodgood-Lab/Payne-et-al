function plotFiringRateMetrics_v1_20250521(inputData, settings)
    % Generates plots related to spatial metrics
    % Written by Anja Payne
    % Last Modified: 05/21/2025
    
    % Inputs:
    %   1) inputData: the matlab structure where the spatial metrics are
    %      saved
    %   2) settings: the settings file where the settings for the spatial
    %      metrics are saved

    % Outputs:
    %   figure 1: trial-by-trial average field size
    %   figure 2: trial-by-trial average field number
    %   figure 3: trial-by-trial avrage information
    %   figure 4: trial-by-trial average sparsity
    
    close all;
    
    if strcmp(settings.rateMaps.plot.display, 'yes') == 1;
        % Get the folders to save plots into
        mainFolder = uigetdir('C:\', 'Please select the folder you would like plots saved into.');
        figureFileInfo = getFigureFolders(mainFolder);
    
        % Have the user select which plots they want to generate
        listOfFigures = getFiguresToPlot();
    
        %% Step 1: Plot the average size, information, and sparsity
        if ismember(1, listOfFigures) == 1; 
            for iGenotype = 1:length(fieldnames(inputData));
                genotypes = fieldnames(inputData); 
                genotypeData = inputData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;
                
                allSize = []; allInfo = []; allSparsity = []; 
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
                                
                                directions = fieldnames(FRdata{iAnimal}(iCluster).rateMaps.rateMap);
                                for iDir = 1:length(directions); 
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);

                                    % Note from AP: outliers were excluded at the suggestion of BB after first
                                    % confirming it did not change statistical outcome
                                    if FRdata{iAnimal}(iCluster).spatialMetrics.PFsize.(directions{iDir}) < 128; 
                                        % Concatenate data
                                        allSize = [allSize, FRdata{iAnimal}(iCluster).spatialMetrics.PFsize.(directions{iDir})];    
                                    end
                                        allInfo = [allInfo, FRdata{iAnimal}(iCluster).spatialMetrics.info.(directions{iDir})];    
                                        allSparsity = [allSparsity, FRdata{iAnimal}(iCluster).spatialMetrics.sparsity.(directions{iDir})]; 
                                    
                                end
                            end
                        end
                    end
                end
                
                size_bothGenotypes{iGenotype} = allSize;
                info_bothGenotypes{iGenotype} = allInfo;
                sparsity_bothGenotypes{iGenotype} = allSparsity;
            end
            
            %%% Plot average PF size
            figures.size = figure(1); clf; hold on;
            plt_WT = cdfplot(size_bothGenotypes{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(size_bothGenotypes{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Average PF Size (cm.)'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(size_bothGenotypes{1})
            KOmean = nanmean(size_bothGenotypes{2})
            WT_SEM = nanstd(size_bothGenotypes{1})/sqrt(length(size_bothGenotypes{1}))
            KO_SEM = nanstd(size_bothGenotypes{2})/sqrt(length(size_bothGenotypes{2}))
            [f, x] = ecdf(size_bothGenotypes{1});
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(size_bothGenotypes{2});
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % If data is not normal, perform kstest
            p = ranksum(size_bothGenotypes{1}, size_bothGenotypes{2});
            pValueDisplay = ['P = ', num2str(p), ' using kstest with WT N = ', ...
                num2str(length(size_bothGenotypes{1})), ' and KO N = ', num2str(length(size_bothGenotypes{2}))];
            display(pValueDisplay); 
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureFileInfo.filePath.population;
            figureSettings.name = 'AverageSize';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureFileInfo.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.size, figureSettings);
            
            %%%%% Plot info %%%%%
            
            figures.info = figure(2); clf; hold on;
            plt_WT = cdfplot(info_bothGenotypes{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(info_bothGenotypes{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Information (bits/spike)'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(info_bothGenotypes{1}); 
            KOmean = nanmean(info_bothGenotypes{2});
            WT_SEM = nanstd(info_bothGenotypes{1})/sqrt(length(info_bothGenotypes{1})); 
            KO_SEM = nanstd(info_bothGenotypes{2})/sqrt(length(info_bothGenotypes{2})); 
            [f, x] = ecdf(info_bothGenotypes{1});
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(info_bothGenotypes{2});
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(info_bothGenotypes{1})
            [h, pUnif_KO] = adtest(info_bothGenotypes{2})
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(info_bothGenotypes{1}, info_bothGenotypes{2});
            pValueDisplay = ['P = ', num2str(p), ' using kstest with WT N = ', ...
                num2str(length(info_bothGenotypes{1})), ' and KO N = ', num2str(length(info_bothGenotypes{2}))];
            display(pValueDisplay);          
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureFileInfo.filePath.population;
            figureSettings.name = 'Information';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureFileInfo.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.info, figureSettings);
            
            %%%%% Plot sparsity %%%%%
            
            figures.sparsity = figure(3); clf; hold on;
            plt_WT = cdfplot(sparsity_bothGenotypes{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(sparsity_bothGenotypes{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Sparsity'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(sparsity_bothGenotypes{1}); 
            KOmean = nanmean(sparsity_bothGenotypes{2});
            WT_SEM = nanstd(sparsity_bothGenotypes{1})/sqrt(length(sparsity_bothGenotypes{1})); 
            KO_SEM = nanstd(sparsity_bothGenotypes{2})/sqrt(length(sparsity_bothGenotypes{2})); 
            [f, x] = ecdf(sparsity_bothGenotypes{1});
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(sparsity_bothGenotypes{2});
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(sparsity_bothGenotypes{1})
            [h, pUnif_KO] = adtest(sparsity_bothGenotypes{2})
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(sparsity_bothGenotypes{1}, sparsity_bothGenotypes{2});
            pValueDisplay = ['P = ', num2str(p), ' using kstest with WT N = ', ...
                num2str(length(sparsity_bothGenotypes{1})), ' and KO N = ', num2str(length(sparsity_bothGenotypes{2}))];
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureFileInfo.filePath.population;
            figureSettings.name = 'Sparsity';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureFileInfo.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.sparsity, figureSettings);
        end
            
        %% Step 2: Trial-by-trial plots
        if ismember(2, listOfFigures) == 1; 
            % Get list of genotypes to plot over
            if strcmp(settings.rateMaps.plot.genotypes, 'all') == 1; 
                genotypeList = 1:length(fieldnames(inputData));
            else
                genotypeList = settings.rateMaps.plot.genotypes; 
            end
            for iGenotype = genotypeList;
                genotypes = fieldnames(inputData); 
                genotypeData = inputData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;
                meanSize = []; meanNumber = []; meanInfo = []; meanSparsity = [];
                
                % Get animals to plot over
                if strcmp(settings.rateMaps.plot.animals, 'all') == 1; 
                    animalList = 1:length(FRdata); 
                else
                    animalList = settings.rateMaps.plot.animals; 
                end
                
                for iAnimal = animalList; 
                    % Skip if empty
                    if isempty(FRdata{iAnimal}) == 1; 
                        continue
                    else
                        
                        % Get cells to plot over
                        if strcmp(settings.rateMaps.plot.cells, 'all') == 1; 
                            [~,n] = size(FRdata{iAnimal});
                            cellList = 1:n; 
                        else
                            cellList = settings.rateMaps.plot.cells; 
                        end
                        
                        for iCluster = cellList; 
                            % Skip if empty
                            if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                                display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                                continue
                            else
                                
                                % Get directions to plot over
                                directions = fieldnames(FRdata{iAnimal}(iCluster).rateMaps.rateMap);
                                if strcmp(settings.rateMaps.plot.direction, 'all') == 1; 
                                    directionList = 1:length(directions); 
                                else
                                    directionList = settings.phasePrecession.plot.direction; 
                                end
                                
                                for iDir = directionList;
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                                     
                                    % Get the mean values
                                    % Note that we only get the mean size
                                    % for the field with the largest FR
                                    meanSize = [meanSize, nanmean(FRdata{iAnimal}(iCluster).spatialMetrics.controls.byTrial.PFsize.(directions{iDir}){1})]; 
                                    meanNumber = [meanNumber, nanmean(FRdata{iAnimal}(iCluster).spatialMetrics.controls.byTrial.PFnumber.(directions{iDir}))]; 
                                    meanInfo = [meanInfo, nanmean(FRdata{iAnimal}(iCluster).spatialMetrics.controls.byTrial.info.(directions{iDir}))]; 
                                    meanSparsity = [meanSparsity, nanmean(FRdata{iAnimal}(iCluster).spatialMetrics.controls.byTrial.sparsity.(directions{iDir}))]; 
                                end
                            end
                        end
                    end
                end
                
                meanSize_bothGenotypes{iGenotype} = meanSize;
                meanNumber_bothGenotypes{iGenotype} = meanNumber;
                meanInfo_bothGenotypes{iGenotype} = meanInfo;
                meanSparsity_bothGenotypes{iGenotype} = meanSparsity;
            end
            
            %%% Plot average field size
            figures.populationMeanSize = figure(1); clf; hold on;
            plt_WT = cdfplot(meanSize_bothGenotypes{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(meanSize_bothGenotypes{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Average Field Size'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(meanSize_bothGenotypes{1}); 
            KOmean = nanmean(meanSize_bothGenotypes{2});
            WT_SEM = nanstd(meanSize_bothGenotypes{1})/sqrt(length(meanSize_bothGenotypes{1})); 
            KO_SEM = nanstd(meanSize_bothGenotypes{2})/sqrt(length(meanSize_bothGenotypes{2})); 
            [f, x] = ecdf(meanSize_bothGenotypes{1});
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(meanSize_bothGenotypes{2});
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(meanSize_bothGenotypes{1})
            [h, pUnif_KO] = adtest(meanSize_bothGenotypes{2})
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(meanSize_bothGenotypes{1}, meanSize_bothGenotypes{2});
            pValueDisplay = ['P = ', num2str(p), ' using kstest with WT N = ', ...
                num2str(length(meanSize_bothGenotypes{1})), ' and KO N = ', num2str(length(meanSize_bothGenotypes{2}))];            
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureFileInfo.filePath.population;
            figureSettings.name = 'TrialByTrial_averageSize';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureFileInfo.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationMeanSize, figureSettings);
            
            %%% Plot average number of fields
            figures.populationNumberOfFields = figure(2); clf; hold on;
            plt_WT = cdfplot(meanNumber_bothGenotypes{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(meanNumber_bothGenotypes{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Average Information'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(meanNumber_bothGenotypes{1}); 
            KOmean = nanmean(meanNumber_bothGenotypes{2});
            WT_SEM = nanstd(meanNumber_bothGenotypes{1})/sqrt(length(meanNumber_bothGenotypes{1})); 
            KO_SEM = nanstd(meanNumber_bothGenotypes{2})/sqrt(length(meanNumber_bothGenotypes{2})); 
            [f, x] = ecdf(meanNumber_bothGenotypes{1});
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(meanNumber_bothGenotypes{2});
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(meanNumber_bothGenotypes{1})
            [h, pUnif_KO] = adtest(meanNumber_bothGenotypes{2})
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(meanNumber_bothGenotypes{1}, meanNumber_bothGenotypes{2});
            pValueDisplay = ['P = ', num2str(p), ' using kstest'];
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureFileInfo.filePath.population;
            figureSettings.name = 'TrialByTrial_averageNumber';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureFileInfo.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationNumberOfFields, figureSettings);
            
            %%% Plot average information
            figures.populationMeanInfo = figure(3); clf; hold on;
            plt_WT = cdfplot(meanInfo_bothGenotypes{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(meanInfo_bothGenotypes{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Average Sparsity'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(meanInfo_bothGenotypes{1}); 
            KOmean = nanmean(meanInfo_bothGenotypes{2});
            WT_SEM = nanstd(meanInfo_bothGenotypes{1})/sqrt(length(meanInfo_bothGenotypes{1})); 
            KO_SEM = nanstd(meanInfo_bothGenotypes{2})/sqrt(length(meanInfo_bothGenotypes{2})); 
            [f, x] = ecdf(meanInfo_bothGenotypes{1});
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(meanInfo_bothGenotypes{2});
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(meanInfo_bothGenotypes{1})
            [h, pUnif_KO] = adtest(meanInfo_bothGenotypes{2})
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(meanInfo_bothGenotypes{1}, meanInfo_bothGenotypes{2});
            pValueDisplay = ['P = ', num2str(p), ' using kstest'];
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureFileInfo.filePath.population;
            figureSettings.name = 'TrialByTrial_averageInfo';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureFileInfo.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationMeanInfo, figureSettings);
            
            %%% Plot average sparsity
            figures.populationMeanSparsity = figure(4); clf; hold on;
            plt_WT = cdfplot(meanSparsity_bothGenotypes{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(meanSparsity_bothGenotypes{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Average Number of Fields'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(meanSparsity_bothGenotypes{1}); 
            KOmean = nanmean(meanSparsity_bothGenotypes{2});
            WT_SEM = nanstd(meanSparsity_bothGenotypes{1})/sqrt(length(meanSparsity_bothGenotypes{1})); 
            KO_SEM = nanstd(meanSparsity_bothGenotypes{2})/sqrt(length(meanSparsity_bothGenotypes{2})); 
            [f, x] = ecdf(meanSparsity_bothGenotypes{1});
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(meanSparsity_bothGenotypes{2});
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(meanSparsity_bothGenotypes{1})
            [h, pUnif_KO] = adtest(meanSparsity_bothGenotypes{2})
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(meanSparsity_bothGenotypes{1}, meanSparsity_bothGenotypes{2});
            pValueDisplay = ['P = ', num2str(p), ' using kstest'];
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureFileInfo.filePath.population;
            figureSettings.name = 'TrialByTrial_averageSparsity';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureFileInfo.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationMeanSparsity, figureSettings);
            
        end
        
    end
end
         
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selectedPlots = getFiguresToPlot()
   % Define the list of available plots
    plotOptions = {'Plot spatial metrics', ...
        'Trial-by-trial spatial metrics',...
        }; 
    
    % Display a dialog box to select the plots
    selectedPlots = listdlg('ListString', plotOptions, ...
    	'SelectionMode', 'multiple', ...
    	'PromptString', 'Select the plots you want to generate:', ...
    	'ListSize', [300, 200]);
end

function figureFileInfo = getFigureFolders(mainFolder)

    % For the population plots
    figureFileInfo.fileNameBase.population = 'spatialMetrics_populationPlots';
    figureFileInfo.filePath.population = getMostRecentFilePath_v1_20240723(figureFileInfo.fileNameBase.population, '', mainFolder);

end

     