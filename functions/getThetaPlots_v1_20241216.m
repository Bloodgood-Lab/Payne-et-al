function data = getThetaPlots_v1_20241216(data)
    % Gets the theta angle of spikes
    % Written by Anja Payne
    % Last Modified: 12/16/2024

    % Inputs:
    %   1) data: the matlab structure where the theta modulation data is
    %      stored
    
    % Outputs:
    %   1) figure 1: a scatter plot on a rose diagram showing the MVL and
    %      preferred direction for all WT and KO cells

    % Steps:
    %   1) Plot figure 1: the population scatter plot
    
    if strcmp(settings.plot, 'yes') == 1; 
        
        % Have the user select which plots they want to generate
        listOfFigures = getFiguresToPlot();
        
        if ismember(1, listOfFigures) == 1; 

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
        
        if ismember(2, listOfFigures) == 1; 

            close all; figure(1); 
            subplot(2,1,1); hold on; 
            lfp_external = thetaData_allTrials.lfp;
            signal_filtered_external = thetaData_allTrials.signal_filtered;
            tSp_external = thetaData_allTrials.tSp; 
            plot(dataForPhaseLocking_allTrials.LFPtimes, lfp_external, 'k')
            plot(dataForPhaseLocking_allTrials.LFPtimes, signal_filtered_external, 'b'); 
            scatter(tSp_external, 0.5*ones(1, length(tSp_external)), '.r');
            subplot(2,1,2); 
            scatter(tSp_external, thetaData_allTrials.allPhs, '.r');
        end
    end
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectedPlots = getFiguresToPlot()
   % Define the list of available plots
    plotOptions = {'preferred phase vs. MVL plotted as vectors on rose plot',...
        'LFP with theta-filtered LFP and spikes'}; 
    
    % Display a dialog box to select the plots
    selectedPlots = listdlg('ListString', plotOptions, ...
    	'SelectionMode', 'multiple', ...
    	'PromptString', 'Select the plots you want to generate:', ...
    	'ListSize', [300, 100]);
end
    
    
    