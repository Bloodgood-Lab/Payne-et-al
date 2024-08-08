function data = getThetaModulation_v1_20240806(data, settings)
    % Gets the theta angle of spikes
    % Written by Anja Payne
    % Last Modified: 08/06/2024

    % Inputs:
    %   1) data: the matlab structure where the in-field spikes-by-trial
    %      are saved
    %   2) settings: the settings file where the frequency band information
    %      for theta is saved
    
    % Outputs:
    %   1) data: the structure with the theta phases appended

    % Steps:
    %   1) Loop through genotypes, animals, and cells and get the theta
    %   phase
    %   2) Ask the user if they want to save the newly generated data
    
    %% Step 1: Get the theta phases
    directoryPath = ''; tetrode = NaN; 
    for iGenotype = 1:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype}); 
        
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
                        newDirectoryPath = FRdata{iAnimal}(iCluster).metaData.directory{1};
                        clusterFileName = FRdata{iAnimal}(iCluster).metaData.fileName{1};
                        newTetrode = clusterFileName(14); 
                        if strcmp(directoryPath, newDirectoryPath) == 1 && tetrode == newTetrode;
                            % If the directory and tetrode are the
                            % same, the CSC for this data has already
                            % been loaded, skip it
                        elseif strcmp(directoryPath, newDirectoryPath) == 0 && tetrode ~= newTetrode;
                            directoryPath = newDirectoryPath; 
                            tetrode = newTetrode; 
                            csc_file = [directoryPath, '\CSC', tetrode, '.ncs'];
                            csc_data = readCSC(csc_file); 
                            csc_timepoints_lowSampling = csc_data.ts/1000; 
                            csc_samples = csc_data.samp; csc_samples = csc_samples(1:16:end); 
                            csc_timepoints_raw = linspace(min(csc_timepoints_lowSampling), max(csc_timepoints_lowSampling), length(csc_samples)); 
                            csc_timepoints = csc_timepoints_raw - min(csc_timepoints_raw); 
                            csc_sampFreq = csc_data.sampFreq(1)/16; 
                        end

                        directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode);
                        for iDir = 1:length(directions);
                            if strcmp(directions(iDir), 'cw') == 1; 
                                spikesByDirection = FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.cw; 
                            elseif strcmp(directions(iDir), 'ccw') == 1; 
                                spikesByDirection = FRdata{iAnimal}(iCluster).inField.inFieldSpkTimes.ccw; 
                            end

                            for iField = 1:length(spikesByDirection); 
                                for iTrial = 1:length(spikesByDirection{iField});
                                    spikesByTrial = spikesByDirection{iField}{iTrial}; 
                                    if isempty(spikesByTrial) == 0; 
                                        % Define input data to theta
                                        % phase locking
                                        dataForPhaseLocking.spikeTimes = round(spikesByTrial);
                                        dataForPhaseLocking.LFPsamples = csc_samples;
                                        dataForPhaseLocking.LFPtimes = round(csc_timepoints);
                                        dataForPhaseLocking.samplingFreq = csc_sampFreq; 
                                        thetaBand = settings.theta.frequencyBand; 
                                        % Get the theta phase locking
                                        thetaData = getPhaseLocking(dataForPhaseLocking, thetaBand);
                                        allTrialPhases{iTrial} = thetaData.allPhs;

                                    elseif isempty(spikesByTrial) == 1; 
                                        allTrialPhases{iTrial} = []; 
                                    end
                                end

                                if strcmp(directions(iDir), 'cw') == 1;
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.phases.cw = allTrialPhases;
                                elseif strcmp(directions(iDir), 'ccw') == 1; 
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.phases.ccw = allTrialPhases;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% Step 2: Save
    saveFile_v1_20240718(data, settings, 'theta') 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function outputData = getPhaseLocking(inputData, band)
        % Gets the theta phase of firing for each spike
        % Code from Sunandha Srikanth in S. Leutgeb lab
        % Inputs: 
        %   1) inputData.spikeTimes: spike times
        %   2) inputData.LFPsamples: LFP samples
        %   3) inputData.LFPtimes: LFP timestamps
        %   4) inputData.samplingFreq: sampling frequency of LFP
        %   5) band: the frequency band to use for extracting the
        %      theta-filtered LFP
        % Output: 
        %   1) outputData.allPhs: the phase for each spike
        %   2) outputData.mvl: preferred phase (mean resultant vector of 
        %      all spike phases) in radians
        %   3) outputData.dir: preferred angle in degrees
        %   4) outputData.p: statistical p-value to check if MVL is
        %      significant
        %   5) outputData.z: z-value obtained from stats
       
        % Define inputs
        tSp = inputData.spikeTimes;
        lfp = inputData.LFPsamples;
        lfp_ts = inputData.LFPtimes;
        Fs = inputData.samplingFreq; 
        
        Wn_theta = [band(1)/(Fs/2) band(2)/(Fs/2)]; % normalized by the nyquist frequency
        [btheta, atheta] = butter(3,Wn_theta);
        signal_filtered = filtfilt(btheta, atheta, lfp);
        thetaPhs = rad2deg(angle(hilbert(signal_filtered)));

        [match_time_spikes, match_time, ~] = intersect(lfp_ts, tSp);

        outputData.allPhs = deg2rad(thetaPhs(match_time));
        outputData.mvl = circ_r(deg2rad(thetaPhs(match_time)));
        outputData.dir = rad2deg(circ_mean(deg2rad(thetaPhs(match_time))));
        [outputData.p, outputData.z] = circ_rtest(deg2rad(thetaPhs(match_time)));
    end

end
    
    
    