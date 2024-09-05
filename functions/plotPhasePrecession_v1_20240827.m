function plotPhasePrecession_v1_20240827(data, settings)
    % Generates plots related to phase precession analysis
    % Written by Anja Payne
    % Last Modified: 08/28/2024

    % Inputs:
    %   1) data: the matlab structure where the in-field theta phases
    %      by-trial are saved
    %   2) settings: the settings file where the settings for the phase
    %      precession calculations are saved

    % Outputs:
    %   1) plots saved in the folder specified by the user
    
    fileNameBase = 'phasePrecession_spikesAndSlopes';
    if strcmp(settings.phasePrecession.plot, 'yes') == 1;
        % Ask the user to select the folder to save figures into
        figureSettings.filePath = getMostRecentFilePath_v1_20240723(fileNameBase, 'Select directory to save figures into');

        for iGenotype = 1:length(fieldnames(data.cellData));
            genotypes = fieldnames(data.cellData);
            genotypeData = data.cellData.(genotypes{iGenotype});

            % Run analysis for high-firing cells only
            FRdata = genotypeData.highFiring;
            for iAnimal = 1%:length(FRdata);
                % Skip if empty
                if isempty(FRdata{iAnimal}) == 1;
                    continue
                else
                    [~,n] = size(FRdata{iAnimal});
                    for iCluster = 2%:n;
                        % Skip if empty
                        if isempty(FRdata{iAnimal}(iCluster).metaData) == 1;
                            display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                            continue
                        else
                            display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                            % Assign variables based on running direction
                            directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode);
                            for iDir = 1:length(directions);
                                outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir));
                                spkPhs = outputData.spkPhsForPlot; spkPos = outputData.spkPosForPlot; binnedSpkPos = outputData.binnedSpkPos; 
                                slopeMedian = outputData.slopeMedian; allSlopes = outputData.allSlopes; lineFit = outputData.lineFit; 
                                
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
                                        figures.allTrialsFig = figure(1); set(figures.allTrialsFig, 'Position', [100, 200, 1800, 800]);
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
                                        ylim([0, 4*pi]); 
                                        
                                        % If there are enough spatial bins
                                        if nanmax(binnedSpkPos{iField}{iTrial})-nanmin(binnedSpkPos{iField}{iTrial}) < settings.phasePrecession.spatialBinThreshold;                                             continue;
                                            continue;
                                        else
                                            % If there were spikes in-field that trial
                                            if nanmax(spkPhs{iField}{iTrial}) > settings.phasePrecession.spikeThreshold; 
                                                % Get the fit values
                                                plotSlope = lineFit{iField}(iTrial).lin.Alpha;
                                                yOffset = lineFit{iField}(iTrial).lin.Phi0;
                                                allPlotSlopes = [allPlotSlopes, plotSlope]; 
                                                allYoffsets = [allYoffsets, yOffset]; 
                                                maxX = nanmax([maxX, nanmax(spkPos{iField}{iTrial})]);
                                                
                                                % Plot that trial's scatter plot
                                                xPlot = [0:max(spkPos{iField}{iTrial})-min(spkPos{iField}{iTrial})+1];
                                                yPlot1 = (xPlot*(plotSlope) + yOffset) - pi;
                                                yPlot2 = (xPlot*(plotSlope) + yOffset) + pi;
                                                plot(xPlot, yPlot1, 'r');
                                                plot(xPlot, yPlot2, 'r');     
                                                if length(xPlot>1); xlim([min(xPlot), max(xPlot)]); end
                                                ylim([0, 4*pi]); 
                                                
                                                % Plot the overlay of all slopes
                                                subplot(subplotNumbRows, 13, allLinePlotRange); hold on;
                                                plot(xPlot, yPlot1, 'k');
                                            end
                                        end
                                    end   
                                    % Clean up the scatter plots
                                    cleanUpPlot(settings, slopeMedian(iField), nanmean(allSlopes{iField}));
                                    
                                    % Clean up the slope overlay plot and
                                    % plot the average slope
                                    inputData.subplotNumbRows = subplotNumbRows;
                                    inputData.allLinePlotRange = allLinePlotRange;
                                    inputData.maxX = maxX; inputData.allPlotSlopes = allPlotSlopes; 
                                    inputData.yOffset = allYoffsets; 
                                    inputData.median = slopeMedian(iField);
                                    inputData.mean = nanmean(allSlopes{iField});
                                    cleanUpAndPlotOverlay(inputData);
                                    
                                    % Plot a histogram of all slopes
                                    plotHistogram(subplotNumbRows, histoPlotRange, allSlopes{iField}, slopeMedian(iField));
                    
                                    % Save the figure
                                    figureSettings.name = [genotypes{iGenotype}, '_Animal', num2str(iAnimal), '_Cluster', ...
                                        num2str(iCluster), '_', directions{iDir}, '_Field', num2str(iField)]; 
                                    figureSettings.appendedFolder.binary = 'yes'; 
                                    figureSettings.appendedFolder.name = fileNameBase;
                                    figureSettings.fileTypes = {'fig', 'tiff'};
                                    saveFigure_v1_20240902(figures.allTrialsFig, figureSettings)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end         
         
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function cleanUpPlot(settings, slopeMedian, slopeMean)
    han = axes('Position', [0.125 0.125 0.8 0.8], 'Visible', 'off');
    han.YLabel.Visible = 'on';
    ylabel(han, 'Theta Phase (degrees)', 'FontSize', 20);
    han = axes('Position', [0.065 0.10 0.9 0.9], 'Visible', 'off');
    han.XLabel.Visible = 'on';
    xlabel(han, 'Linear Position', 'FontSize', 20);
    annotation('textbox', [0.02, 0, 0.2, 0.03], 'String', ...
        ['Data File: ', settings.dataSavePath], 'FitBoxToText', 'on', 'BackgroundColor', 'none', ...
        'EdgeColor', 'none');
    annotation('textbox', [0.75, 0.95, 0.2, 0.03], 'String', ...
        ['Median Slope: ', num2str(slopeMedian)], 'FitBoxToText', 'on', ...
        'BackgroundColor', 'none', 'EdgeColor', 'none');
    annotation('textbox', [0.87, 0.95, 0.2, 0.03], 'String', ...
        ['Mean Slope: ', num2str(slopeMean)], 'FitBoxToText', 'on', ...
        'BackgroundColor', 'none', 'EdgeColor', 'none');
    set(gcf, 'PaperUnits', 'Inches', 'PaperPositionMode', 'auto');
end

function cleanUpAndPlotOverlay(inputData)
    subplot(inputData.subplotNumbRows, 13, inputData.allLinePlotRange); hold on;
    xAverageLine = 0:inputData.maxX;
    inputData.allPlotSlopes(isinf(inputData.allPlotSlopes)==1) = NaN;
    inputData.yOffset(isinf(inputData.yOffset)==1) = NaN;
    yMedianLine = xAverageLine*inputData.median + nanmean(inputData.yOffset)-pi;
    yMeanLine = xAverageLine*inputData.mean + nanmean(inputData.yOffset)-pi;
    medianLine = plot(xAverageLine, yMedianLine, 'r', 'LineWidth', 2);
    meanLine = plot(xAverageLine, yMeanLine, 'b', 'LineWidth', 2);
    legend([medianLine, meanLine], {'Median', 'Mean'}, 'Location', 'northwest');
    legend('boxoff');
    ylim([0, 4*pi]); 
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
                    
                    
                    