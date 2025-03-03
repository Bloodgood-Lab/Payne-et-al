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
    %   1) figure 1: the LFP, theta-filtered LFP, and spikes
    %   2) figure 2: the rose plots showing phase preference of all spikes
    %      plotted for each individual neuron
    
    %   1) figure 1: a scatter plot on a rose diagram showing the MVL and
    %      preferred direction for all WT and KO cells

    % Steps/Figures:
    %   1) For indicated animals/cells, get the LFP, theta-filtered LFP,
    %      and spikes in one plot. In another, plot the spike phases. 
    
    %   1) Plot figure 1: the population scatter plot
    
    
    
    % Get the folders to save plots into
    mainFolder = uigetdir('C:\', 'Please select the folder you would like phase precession plots saved into.');
    figureSettings = getFigureFolders(mainFolder); 
    close all; 
    
    if strcmp(settings.plot.display, 'yes') == 1; 
        
        % Have the user select which plots they want to generate
        listOfFigures = getFiguresToPlot();
        
        %% Figure 1: LFP with spikes
        if ismember(1, listOfFigures) == 1; 
            % Get list of genotypes to plot over
            if strcmp(settings.plot.genotypes, 'all') == 1; 
                genotypeList = 1:length(fieldnames(data.cellData));
            else
                genotypeList = settings.plot.genotypes; 
            end
            for iGenotype = genotypeList;
                genotypes = fieldnames(data.cellData); 
                genotypeData = data.cellData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;
                
                % Get animals to plot over
                if strcmp(settings.plot.animals, 'all') == 1; 
                    animalList = 1:length(FRdata); 
                else
                    animalList = settings.plot.animals; 
                end
                
                for iAnimal = animalList; 
                    % Skip if empty
                    if isempty(FRdata{iAnimal}) == 1; 
                        continue
                    else
                        
                        % Get cells to plot over
                        if strcmp(settings.plot.cells, 'all') == 1; 
                            [~,n] = size(FRdata{iAnimal});
                            cellList = 1:n; 
                        else
                            cellList = settings.plot.cells; 
                        end
                        
                        for iCluster = cellList; 
                            % Skip if empty
                            if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                                display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                                continue
                            else
                                
                                % Get directions to plot over
                                directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original);
                                if strcmp(settings.plot.direction, 'all') == 1; 
                                    directionList = 1:length(directions); 
                                else
                                    directionList = settings.plot.direction; 
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

                                        %subplot(2,1,2); 
                                        %scatter(tSp_external, thetaData_allTrials.allPhs, '.r');

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
        
        %% Figure 2: Individual rose plots

        if ismember(2, listOfFigures) == 1; 
            % Get list of genotypes to plot over
            if strcmp(settings.plot.genotypes, 'all') == 1; 
                genotypeList = 1:length(fieldnames(data.cellData));
            else
                genotypeList = settings.plot.genotypes; 
            end
            for iGenotype = genotypeList;
                genotypes = fieldnames(data.cellData); 
                genotypeData = data.cellData.(genotypes{iGenotype}); 
                FRdata = genotypeData.highFiring;

                % Get animals to plot over
                if strcmp(settings.plot.animals, 'all') == 1; 
                    animalList = 1:length(FRdata); 
                else
                    animalList = settings.plot.animals; 
                end

                for iAnimal = animalList; 
                    % Skip if empty
                    if isempty(FRdata{iAnimal}) == 1; 
                        continue
                    else

                        % Get cells to plot over
                        if strcmp(settings.plot.cells, 'all') == 1; 
                            [~,n] = size(FRdata{iAnimal});
                            cellList = 1:n; 
                        else
                            cellList = settings.plot.cells; 
                        end

                        for iCluster = cellList; 
                            % Skip if empty
                            if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                                display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                                continue
                            else

                                % Get directions to plot over
                                directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original);
                                if strcmp(settings.plot.direction, 'all') == 1; 
                                    directionList = 1:length(directions); 
                                else
                                    directionList = settings.plot.direction; 
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
        
        %% Figure 3: Rose plot scatter plot for all neurons
        
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
        
        %% Figure 4: Preferred Phase
        
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
            scatter(x_WT, y_WT, 50, 'ok'); 
            scatter(x_KO, y_KO, 50, 'og'); 

            % Add the unit circle for reference
            theta = linspace(0, 2*pi, 100);
            plot(cos(theta), sin(theta), 'k-'); % Unit circle in dashed black
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
            
            scatter(xStats_WT, yStats_WT, 200, 'filled', 'k'); 
            scatter(xStats_KO, yStats_KO, 200, 'filled', 'g'); 

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
        'preferred phase for WT and KO populations'}; 
    
    % Display a dialog box to select the plots
    selectedPlots = listdlg('ListString', plotOptions, ...
    	'SelectionMode', 'multiple', ...
    	'PromptString', 'Select the plots you want to generate:', ...
    	'ListSize', [300, 100]);
end
    
    
    