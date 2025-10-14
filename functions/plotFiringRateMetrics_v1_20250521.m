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
    %   figure 1: distribution of mean spatial firing rates
    
    close all;
    
    if strcmp(settings.firingRates.plot.display, 'yes') == 1;
        % Get the folders to save plots into
        mainFolder = uigetdir('C:\', 'Please select the folder you would like phase precession plots saved into.');
        figureFileInfo = getFigureFolders(mainFolder);
    
        % Have the user select which plots they want to generate
        listOfFigures = getFiguresToPlot();
    
        %% Figure 1: Mean and Max Spatial FRs
        if ismember(1, listOfFigures) == 1; 
            for iGenotype = 1:length(fieldnames(inputData));
                genotypes = fieldnames(inputData); 
                genotypeData = inputData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;
            
                meanFR = []; maxFR = [];
                for iAnimal = 1:length(FRdata); 
                    % Skip if empty
                    if isempty(FRdata{iAnimal}) == 1; 
                        continue
                    else
                  
                        for iCluster = 1:length(FRdata{iAnimal});
                            % Skip if empty
                            if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                                display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                                continue
                            else
                                
                                % Get directions to plot over
                                directions = fieldnames(FRdata{iAnimal}(iCluster).rateMaps.rateMap);
                                
                                for iDir = 1:length(directions);
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                                     
                                    % Get the mean and max FRs
                                    meanFR = [meanFR, FRdata{iAnimal}(iCluster).firingRates.mean.(directions{iDir})]; 
                                    maxFR = [maxFR, FRdata{iAnimal}(iCluster).firingRates.max.(directions{iDir})]; 
                                end
                            end
                        end
                    end
                end
                
                meanFR_bothGenotypes{iGenotype} = meanFR;
                maxFR_bothGenotypes{iGenotype} = maxFR;
            end
            
            %%% Plot average FR
            figures.populationMeanFR = figure(1); clf; hold on;
            histogram(meanFR_bothGenotypes{1}, 'BinEdges', [0:0.1:5], 'Normalization',...
                'probability', 'FaceColor', [0.5, 0.5, 0.5]); 
            histogram(meanFR_bothGenotypes{2}, 'BinEdges', [0:0.1:5], 'Normalization',...
                'probability', 'FaceColor', [0.0, 0.5, 0.0]); 
            xlabel('Mean FR (Hz.)'); ylabel('Proportion'); 
            set(gca, 'FontSize', 14); 

            % Save the figure
            figureSettings.filePath = figureFileInfo.filePath.population;
            figureSettings.name = 'meanFR_histogram';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureFileInfo.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationMeanFR, figureSettings);
            
            %%% Plot max FR
            figures.populationMaxFR = figure(2); clf; hold on;
            histogram(maxFR_bothGenotypes{1}, 'BinEdges', [0:1:30], 'Normalization',...
                'probability', 'FaceColor', [0.5, 0.5, 0.5]); 
            histogram(maxFR_bothGenotypes{2}, 'BinEdges', [0:1:30], 'Normalization',...
                'probability', 'FaceColor', [0.0, 0.5, 0.0]); 
            xlabel('Max FR (Hz.)'); ylabel('Proportion'); 
            set(gca, 'FontSize', 14); 

            % Save the figure
            figureSettings.filePath = figureFileInfo.filePath.population;
            figureSettings.name = 'maxFR_histogram';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureFileInfo.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationMaxFR, figureSettings);
        end
        
        %% Figure 2: Get burst index for low- and high-firing cells
        if ismember(2, listOfFigures) == 1; 
            for iGenotype = 1:length(fieldnames(inputData));
                genotypes = fieldnames(inputData); 
                genotypeData = inputData.(genotypes{iGenotype}); 
                highFR = genotypeData.highFiring;
                lowFR = genotypeData.lowFiring;
            
                burstIndices_highFiring = []; burstIndices_lowFiring = []; 
                for iAnimal = 1:length(highFR); 
                    % Skip if  high FR is empty
                    if isempty(highFR{iAnimal}) == 1; 
                        continue
                    else
                  
                        for iCluster = 1:length(highFR{iAnimal});
                            % Skip if empty
                            if isempty(highFR{iAnimal}(iCluster).metaData) == 1; 
                                display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                                continue
                            else
                                
                                % Get directions to plot over
                                directions = fieldnames(highFR{iAnimal}(iCluster).rateMaps.rateMap);
                                
                                for iDir = 1:length(directions);
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                                     
                                    % Concatenate all bursts 
                                    allBursts = []; allSpikes = []; 
                                    for iTrial = 1:length(highFR{iAnimal}(iCluster).spikeTiming.bursts.(directions{iDir})); 
                                       allBursts = [allBursts, highFR{iAnimal}(iCluster).spikeTiming.bursts.(directions{iDir}){iTrial}]; 
                                       allSpikes = [allSpikes; highFR{iAnimal}(iCluster).binnedSpikesByTrial.highVelocityData.binnedSpkTimes.(directions{iDir}){iTrial}]; 
                                    end
                                    
                                    % Remove NaNs
                                    allBursts_withoutNaN = allBursts(~isnan(allBursts));
                                    allSpikes_withoutNaN = allSpikes(~isnan(allSpikes)); 
                                    
                                    % Get burst index
                                    burstIndex = length(allBursts_withoutNaN) / length(allSpikes_withoutNaN); 
                                    burstIndices_highFiring = [burstIndices_highFiring, burstIndex]; 
                                end
                            end
                        end
                    end
                    
                    % Skip if  low FR is empty
                    if isempty(lowFR{iAnimal}) == 1; 
                        continue
                    else
                  
                        for iCluster = 1:length(lowFR{iAnimal});
                            % Skip if empty
                            if isempty(lowFR{iAnimal}(iCluster).metaData) == 1; 
                                display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                                continue
                            else
                                
                                % Get directions to plot over
                                directions = fieldnames(lowFR{iAnimal}(iCluster).rateMaps.rateMap);
                                
                                for iDir = 1:length(directions);
                                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                                     
                                    % Concatenate all bursts 
                                    allBursts = []; allSpikes = []; 
                                    for iTrial = 1:length(lowFR{iAnimal}(iCluster).spikeTiming.bursts.(directions{iDir})); 
                                       allBursts = [allBursts, lowFR{iAnimal}(iCluster).spikeTiming.bursts.(directions{iDir}){iTrial}]; 
                                       allSpikes = [allSpikes; lowFR{iAnimal}(iCluster).binnedSpikesByTrial.highVelocityData.binnedSpkTimes.(directions{iDir}){iTrial}]; 
                                    end
                                    
                                    % Remove NaNs
                                    allBursts_withoutNaN = allBursts(~isnan(allBursts));
                                    allSpikes_withoutNaN = allSpikes(~isnan(allSpikes)); 
                                    
                                    % Get burst index
                                    burstIndex = length(allBursts_withoutNaN) / length(allSpikes_withoutNaN); 
                                    burstIndices_lowFiring = [burstIndices_lowFiring, burstIndex]; 
                                end
                            end
                        end
                    end
                end
                
                burstIndices_highFiring_bothGenotypes{iGenotype} = burstIndices_highFiring;
                burstIndices_lowFiring_bothGenotypes{iGenotype} = burstIndices_lowFiring;
            end
            
            %%% Plot high firing burst index
            figures.populationburstHighFR = figure(3); clf; hold on;
            plt_WT = cdfplot(burstIndices_highFiring_bothGenotypes{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(burstIndices_highFiring_bothGenotypes{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Proportion of Spikes in Burst'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(burstIndices_highFiring_bothGenotypes{1}); 
            KOmean = nanmean(burstIndices_highFiring_bothGenotypes{2});
            WT_SEM = nanstd(burstIndices_highFiring_bothGenotypes{1})/sqrt(length(burstIndices_highFiring_bothGenotypes{1})); 
            KO_SEM = nanstd(burstIndices_highFiring_bothGenotypes{2})/sqrt(length(burstIndices_highFiring_bothGenotypes{2})); 
            [f, x] = ecdf(burstIndices_highFiring_bothGenotypes{1});
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(burstIndices_highFiring_bothGenotypes{2});
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(burstIndices_highFiring_bothGenotypes{1})
            [h, pUnif_KO] = adtest(burstIndices_highFiring_bothGenotypes{2})
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(burstIndices_highFiring_bothGenotypes{1}, burstIndices_highFiring_bothGenotypes{2});
            pValueDisplay = ['P = ', num2str(p), ' using kstest with WT N = ', ...
                num2str(length(burstIndices_highFiring_bothGenotypes{1})), ' and KO N = ', num2str(length(burstIndices_highFiring_bothGenotypes{2}))];            
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureFileInfo.filePath.population;
            figureSettings.name = 'BurstIndex_HighFR';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureFileInfo.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationburstHighFR, figureSettings);
            
            %%% Plot low firing burst index
            figures.populationburstLowFR = figure(4); clf; hold on;
            plt_WT = cdfplot(burstIndices_lowFiring_bothGenotypes{1}); 
            set(plt_WT, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
            plt_KO = cdfplot(burstIndices_lowFiring_bothGenotypes{2}); 
            set(plt_KO, 'Color', [0.0 0.5 0.0], 'LineWidth', 1.5); 
            grid off; 
            xlabel('Proportion of Spikes in Burst'); 
            set(gca, 'FontSize', 14); 

            % Add the mean +/- SEM 
            WTmean = nanmean(burstIndices_lowFiring_bothGenotypes{1}); 
            KOmean = nanmean(burstIndices_lowFiring_bothGenotypes{2});
            WT_SEM = nanstd(burstIndices_lowFiring_bothGenotypes{1})/sqrt(length(burstIndices_lowFiring_bothGenotypes{1})); 
            KO_SEM = nanstd(burstIndices_lowFiring_bothGenotypes{2})/sqrt(length(burstIndices_lowFiring_bothGenotypes{2})); 
            [f, x] = ecdf(burstIndices_lowFiring_bothGenotypes{1});
            [~, idx] = min(abs(x-WTmean)); WT_y = f(idx); 
            [f, x] = ecdf(burstIndices_lowFiring_bothGenotypes{2});
            [~, idx] = min(abs(x-KOmean)); KO_y = f(idx);
            plot([WTmean - WT_SEM, WTmean + WT_SEM], [WT_y, WT_y], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2); 
            plot(WTmean, WT_y, 'o', 'MarkerFaceColor', [0.2, 0.2, 0.2], 'MarkerSize', 6);
            plot([KOmean - KO_SEM, KOmean + KO_SEM], [KO_y, KO_y], 'Color', [0.0, 0.2, 0.0], 'LineWidth', 2); 
            plot(KOmean, KO_y, 'o', 'MarkerFaceColor', [0.0, 0.2, 0.0], 'MarkerSize', 6);

            % Display statistical significance
            % First, check for normality 
            [h, pUnif_WT] = adtest(burstIndices_lowFiring_bothGenotypes{1})
            [h, pUnif_KO] = adtest(burstIndices_lowFiring_bothGenotypes{2})
            if pUnif_WT < 0.05 && pUnif_KO < 0.05
                display('Data is not normal'); 
            end
            % If data is not normal, perform kstest
            [~, p] = kstest2(burstIndices_lowFiring_bothGenotypes{1}, burstIndices_lowFiring_bothGenotypes{2});
            pValueDisplay = ['P = ', num2str(p), ' using kstest with WT N = ', ...
                num2str(length(burstIndices_lowFiring_bothGenotypes{1})), ' and KO N = ', num2str(length(burstIndices_lowFiring_bothGenotypes{2}))]; 
            display(pValueDisplay);
            annotation('textbox', [0.55, 0.01, 0.5, 0.05], 'String', pValueDisplay, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

            % Save the figure
            figureSettings.filePath = figureFileInfo.filePath.population;
            figureSettings.name = 'BurstIndex_LowFR';
            figureSettings.appendedFolder.binary = 'no'; 
            figureSettings.appendedFolder.name = figureFileInfo.fileNameBase.population;
            figureSettings.fileTypes = {'fig', 'tiff'};
            saveFigure_v1_20240902(figures.populationburstLowFR, figureSettings);
        end
        
    end
end
         
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selectedPlots = getFiguresToPlot()
   % Define the list of available plots
    plotOptions = {'Distribution of Mean and Max Spatial FRs',...
        'Burst Index for Low- and High-Firing Cells'}; 
    
    % Display a dialog box to select the plots
    selectedPlots = listdlg('ListString', plotOptions, ...
    	'SelectionMode', 'multiple', ...
    	'PromptString', 'Select the plots you want to generate:', ...
    	'ListSize', [300, 200]);
end

function figureFileInfo = getFigureFolders(mainFolder)

    % For the population plots
    figureFileInfo.fileNameBase.population = 'firingRates_populationPlots';
    figureFileInfo.filePath.population = getMostRecentFilePath_v1_20240723(figureFileInfo.fileNameBase.population, '', mainFolder);

end

     