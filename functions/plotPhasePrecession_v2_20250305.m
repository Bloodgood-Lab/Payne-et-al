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
        % Assign the preferredPhase to WT and KO
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
            %medianSlopes_WT_noNaN = medianSlopes_WT; medianSlopes_WT_noNaN(isnan(medianSlopes_WT_noNaN)) = []; 
            %medianSlopes_KO_noNaN = medianSlopes_KO; medianSlopes_KO_noNaN(isnan(medianSlopes_KO_noNaN)) = []; 
            %WT_y = interp1(sort(medianSlopes_WT_noNaN), linspace(0, 1, length(medianSlopes_WT_noNaN)), WTmean, 'linear', 'extrap');
            %KO_y = interp1(sort(medianSlopes_KO_noNaN), linspace(0, 1, length(medianSlopes_KO_noNaN)), KOmean, 'linear', 'extrap');
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
    end
    
end   
         
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selectedPlots = getFiguresToPlot()
   % Define the list of available plots
    plotOptions = {'LFP with spikes and theta phase',...
        'spikes with slopes and histogram of slopes',...
        'CDF of median slopes for WT and KO populations'}; 
    
    % Display a dialog box to select the plots
    selectedPlots = listdlg('ListString', plotOptions, ...
    	'SelectionMode', 'multiple', ...
    	'PromptString', 'Select the plots you want to generate:', ...
    	'ListSize', [300, 100]);
end

function figureSettings = getFigureFolders(mainFolder)

    % For the LFP with spikes and theta phases
    figureSettings.fileNameBase.LFP = 'phasePrecession_thetaLFPwithSpikesAndSpikePhases';
    figureSettings.filePath.LFP = getMostRecentFilePath_v1_20240723(figureSettings.fileNameBase.LFP, '', mainFolder);
    
    % For the spikes and slopes plots
    figureSettings.fileNameBase.spikesAndSlopes = 'phasePrecession_spikesAndSlopes';
    figureSettings.filePath.spikesAndSlopes = getMostRecentFilePath_v1_20240723(figureSettings.fileNameBase.spikesAndSlopes, '', mainFolder);

    % For the population CDF median slope
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


                    