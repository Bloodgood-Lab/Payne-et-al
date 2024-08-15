function data = getSpikeTimes_v1_20240725(data, settings, processedDataPath)
    % Gets the spike times and appends them to the data structure
    % Written by Anja Payne
    % Last Modified: 07/25/2024

    % Inputs:
    %   1) data: the matlab structure where the file paths to the spike
    %      time files from MClust are stored
    %   2) settings: the matlab settings file where the thresholds are
    %      stored
    
    % Outputs:
    %   1) data: the structure with the spike times in milliseconds
    %      appended

    % Steps:
    %   1) Loop through genotypes, animals, and cells to get the file path,
    %      load the data, and append it to the data structure
    %   2) Ask the user if they want to save the newly generated data
    
    %% Step 1: Append spike times to data structure
    directoryPath = '';
    for iGenotype = 1:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype}); 
        for iAnimal = 1:length(genotypeData); 
            if isempty(genotypeData{iAnimal}) == 1; 
                continue
            else
                [~,n] = size(genotypeData{iAnimal});
                for iCluster = 1:n; 
                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);

                    newDirectoryPath = genotypeData{iAnimal}(iCluster).directory{1};
                    if strcmp(directoryPath, newDirectoryPath) == 1;
                        % Do nothing
                    elseif strcmp(directoryPath, newDirectoryPath) == 0;
                        directoryPath = newDirectoryPath; 
                        % Load the CSC data
                        CSC_file = [directoryPath, '\CSC1.ncs']; 
                        csc_data = readCSC(CSC_file); 
                        csc_timepoints_lowSampling = csc_data.ts; 
                        csc_samples = csc_data.samp; 
                        csc_samples_points = linspace(1, length(csc_timepoints_lowSampling), length(csc_samples)); 
                        csc_timepoints_lowSampling_points = [1:length(csc_timepoints_lowSampling)]; 
                        csc_timepoints = interp1(csc_timepoints_lowSampling_points, csc_timepoints_lowSampling, csc_samples_points, 'linear'); 
                        csc_timepoints_msec = csc_timepoints/1000;
                        ephys_t0 = csc_timepoints_msec(1);
                    end
                    % Load the spike data for the relevant cluster
                    spikeDataPath = strcat(directoryPath, '\', genotypeData{iAnimal}(iCluster).fileName{1}, '.mat');
                    load(spikeDataPath);

                    % Get the spike times in msec and aligned to the beginning 
                    % of the recording and save into the data structure
                    spike_time_msec = TS/10 - ephys_t0;
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).spikeTimes = spike_time_msec; 

                end
            end
        end
    end
    
    %% Step 2: Step 2: Save
    saveFile_v1_20240718(processedDataPath, data, settings, 'spikeTimes') 
