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
                                if strcmp(directions(iDir), 'cw') == 1;
                                    spkPhs = FRdata{iAnimal}(iCluster).phasePrecession.phsInput.cw;
                                    spkPos = FRdata{iAnimal}(iCluster).phasePrecession.posInput.cw;
                                    lineFit = FRdata{iAnimal}(iCluster).phasePrecession.fitInfo.cw;
                                    binnedSpkPos = FRdata{iAnimal}(iCluster).inField.inFieldBinnedSpkPos.cw;
                                    slopeMedian = FRdata{iAnimal}(iCluster).phasePrecession.medianSlope.cw;
                                elseif strcmp(directions(iDir), 'ccw') == 1;
                                    spkPhs = FRdata{iAnimal}(iCluster).phasePrecession.phsInput.ccw;
                                    spkPos = FRdata{iAnimal}(iCluster).phasePrecession.posInput.ccw;
                                    lineFit = FRdata{iAnimal}(iCluster).phasePrecession.fitInfo.ccw;
                                    binnedSpkPos = FRdata{iAnimal}(iCluster).inField.inFieldBinnedSpkPos.ccw;
                                    slopeMedian = FRdata{iAnimal}(iCluster).phasePrecession.medianSlope.ccw;
                                end
                                % Loop through fields
                                if strcmp(settings.phasePrecession.fieldsToAnalyze, 'all fields') == 1;
                                    numFieldsToAnalyze = length(spkPhs);
                                elseif strcmp(settings.phasePrecession.fieldsToAnalyze, 'best field') == 1;
                                    numFieldsToAnalyze = 1;
                                end
                                for iField = 1:numFieldsToAnalyze;
                                    % Loop through all the trials
                                    subplotCount = 1; close all; 
                                    for iTrial = 1:length(spkPhs{iField});
                                        figures.allTrialsFig = figure(1); set(figures.allTrialsFig, 'Position', [100, 200, 1800, 800]);
                                        if iTrial == 1; clf(figures.allTrialsFig); else hold on; end;
                                        subplot(ceil(length(spkPhs{iField})/10),10,subplotCount); hold on;
                                        scatter(spkPos{iField}{iTrial}, spkPhs{iField}{iTrial}, 200, '.k');
                                        set(gca, 'FontSize', 12);
                                        title(num2str(iTrial));
                                        % If there are enough spatial bins
                                        if nanmax(binnedSpkPos{iField}{iTrial})-nanmin(binnedSpkPos{iField}{iTrial}) < settings.phasePrecession.spatialBinThreshold;                                             continue;
                                            continue;
                                        else
                                            % If there were spikes in-field that trial
                                            if nanmax(spkPhs{iField}{iTrial}) > settings.phasePrecession.spikeThreshold; 
                                                xPlot = [0:max(spkPos{iField}{iTrial})-min(spkPos{iField}{iTrial})+1];
                                                yPlot1 = (xPlot*(lineFit{iField}(iTrial).lin.Alpha)+lineFit{iField}(iTrial).lin.Phi0)-pi;
                                                yPlot2 = (xPlot*(lineFit{iField}(iTrial).lin.Alpha)+lineFit{iField}(iTrial).lin.Phi0)+pi;
                                                plot(xPlot, yPlot1, 'r');
                                                plot(xPlot, yPlot2, 'r');     
                                                if length(xPlot>1); xlim([min(xPlot), max(xPlot)]); end
                                            end
                                        end
                                        subplotCount = subplotCount + 1;
                                    end      
                                    han = axes('Position', [0.125 0.125 0.8 0.8], 'Visible', 'off');
                                    han.YLabel.Visible = 'on';
                                    ylabel(han, 'Theta Phase (degrees)', 'FontSize', 20);
                                    han = axes('Position', [0.065 0.10 0.9 0.9], 'Visible', 'off');
                                    han.XLabel.Visible = 'on';
                                    xlabel(han, 'Normalized Linear Position', 'FontSize', 20);
                                    annotation('textbox', [0.02, 0, 0.2, 0.03], 'String', ...
                                        ['Data File: ', settings.dataSavePath], 'FitBoxToText', 'on', 'BackgroundColor', 'none', ...
                                        'EdgeColor', 'none');
                                    annotation('textbox', [0.9, 0, 0.2, 0.03], 'String', ...
                                        ['Median Slope: ', num2str(slopeMedian)], 'FitBoxToText', 'on', ...
                                        'BackgroundColor', 'none', 'EdgeColor', 'none');
                                    set(gcf, 'PaperUnits', 'Inches', 'PaperPositionMode', 'auto');
                    
                                    % Save the figure
                                    figureSettings.name = [genotypes{iGenotype}, ...
                                        '_Animal', num2str(iAnimal), '_Cluster', ...
                                        num2str(iCluster), '_', directions{iDir}, ...
                                        '_Field', num2str(iField)]; 
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
          %{
            % If there are figures, have the user select the save path and
            % save them there
            figs = findobj('Type', 'figure');
            if isempty(figs) == 0; % if there are figures
                for iFig = 1:length(figs)
                    display('Saving Figures')
                    filePath = getMostRecentFilePath_v1_20240723(fileNameBase);
                    savePathName = filePath{1};
                    saveVersion = str2double(filePath{2}) + 1;
                    saveVersion = sprintf('%02d', saveVersion);
                    saveFile = [fileNameBase, '_v', saveVersion, '_', date];
                    if ~exist(savePathName, 'dir')
                        % Create the folder if it does not exist
                        mkdir(savePathName);
                    end
                    saveas(figs(iFig), [savePathName, '\', saveFile]);
                end
            end
    %}          
                    
                    
                    
                    
                    