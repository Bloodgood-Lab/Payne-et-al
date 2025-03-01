function data = getThetaPlots_v1_20241216(data, settings)
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
                                        figures.LFP = figure(1); hold on;
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
                                        genotypes{iGenotype}
                                        directions(iDir)
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
            
        if ismember(2, listOfFigures) == 1; 

            %% Step 1: Get the population scatter plot
            preferred_phase_WT_rad = deg2rad(data.populationData(1).preferredDirection); % Random phases for WT (0 to 2*pi)
            preferred_phase_KO_rad = deg2rad(data.populationData(2).preferredDirection); % Random phases for KO (0 to 2*pi)
            mean_vector_length_WT = data.populationData(1).MVL  ;
            mean_vector_length_KO = data.populationData(2).MVL  ;

            figure;
            % Plot WT data
            polar(0, 0.8, 'ko');
            hold on;
            for i = 1:length(preferred_phase_WT_rad)
                % Use polar function to plot each point for WT
                polar(preferred_phase_WT_rad(i), mean_vector_length_WT(i), 'ko');
            end

            % Plot KO data
            figure;
            polar(0, 0.8, 'ko');
            hold on;
            for i = 1:length(preferred_phase_KO_rad)
                % Use polar function to plot each point for KO
                polar(preferred_phase_KO_rad(i), mean_vector_length_KO(i), 'go');
            end
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

end

function selectedPlots = getFiguresToPlot()
   % Define the list of available plots
    plotOptions = {'LFP with theta-filtered LFP and spikes', ...
        'preferred phase vs. MVL plotted as vectors on rose plot'}; 
    
    % Display a dialog box to select the plots
    selectedPlots = listdlg('ListString', plotOptions, ...
    	'SelectionMode', 'multiple', ...
    	'PromptString', 'Select the plots you want to generate:', ...
    	'ListSize', [300, 100]);
end
    
    
    