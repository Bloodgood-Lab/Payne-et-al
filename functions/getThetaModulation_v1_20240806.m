function data = getThetaModulation_v1_20240806(data, settings, processedDataPath)
    % Gets the theta angle of spikes
    % Written by Anja Payne
    % Last Modified: 12/04/2024

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
    
    
    data.cellData = data; data = rmfield(data, 'WT'); data = rmfield(data, 'KO');
    %% Step 1: Get the theta phases
    directoryPath = ''; tetrode = NaN; 
    for iGenotype = 1:length(fieldnames(data.cellData));
        genotypes = fieldnames(data.cellData); 
        genotypeData = data.cellData.(genotypes{iGenotype}); 
        populationMVL = []; populationDir = []; populationHistogram = [];
        populationBurstsMVL = []; populationBurstsDir = []; populationBurstsHistogram = [];
        populationSinglesMVL = []; populationSinglesDir = []; populationSinglesHistogram = [];
        
        % Run analysis for high-firing cells only
        FRdata = genotypeData.highFiring;
        for iAnimal = 1:length(FRdata); 
            % Skip if empty
            if isempty(FRdata{iAnimal}) == 1; 
                continue
            else
                [~,n] = size(FRdata{iAnimal});
                for iCluster = 1%:n;
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
                        elseif strcmp(directoryPath, newDirectoryPath) == 0 || tetrode ~= newTetrode;
                            directoryPath = newDirectoryPath; 
                            tetrode = newTetrode; 
                            csc_file = [directoryPath, '\CSC', tetrode, '.ncs'];
                            csc_data = readCSC(csc_file); 
                            csc_timepoints_lowSampling = csc_data.ts/1000; 
                            csc_samples = csc_data.samp; csc_samples = csc_samples(1:16:end); 
                            csc_timepoints_raw = linspace(min(csc_timepoints_lowSampling), max(csc_timepoints_lowSampling), length(csc_samples)); 
                            csc_timepoints = csc_timepoints_raw - min(csc_timepoints_raw); 
                            csc_sampFreq = csc_data.sampFreq(1)/16; 
                            % Define these as inputs for the theta
                            % modulation function and clear variables
                            dataForPhaseLocking.LFPsamples = csc_samples;
                            dataForPhaseLocking.LFPtimes = round(csc_timepoints);
                            dataForPhaseLocking.samplingFreq = csc_sampFreq; 
                            clear csc_data csc_timepoints_lowSampling csc_timepoints_raw csc_samples csc_timepoints
                        end

                        directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original);
                        for iDir = 1:length(directions);
                            outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir), 'theta');
                            spikesByDirection = outputData.spikesByDirection; 
                            inFieldBursts = outputData.inFieldBursts;
                            inFieldSingles = outputData.inFieldSingles;

                            % Loop through fields
                            numFieldsToAnalyze = whichField(settings.theta.fieldsToAnalyze, spikesByDirection);
                            allSpikes_allTrials = []; allBursts_allTrials = []; allSingles_allTrials = []; 
                            allTrialPhases = {}; allTrialPhases_bursts = {}; allTrialPhases_singles = {}; 
                            for iField = 1:numFieldsToAnalyze;
                                for iTrial = 1:length(spikesByDirection{iField});
                                    spikesByTrial = spikesByDirection{iField}{iTrial};
                                    burstsByTrial = inFieldBursts{iField}{iTrial};
                                    singlesByTrial = inFieldSingles{iField}{iTrial};
                                    allSpikes_allTrials = [allSpikes_allTrials; spikesByTrial]; 
                                    allBursts_allTrials = [allBursts_allTrials, burstsByTrial]; 
                                    allSingles_allTrials = [allSingles_allTrials, singlesByTrial]; 
                                    
                                    % Get theta modulation for all spikes
                                    if isempty(spikesByTrial) == 0; 
                                        % Define the rest of the inputs for
                                        % the theta modulation script
                                        dataForPhaseLocking.spikeTimes = round(spikesByTrial);
                                        thetaBand = settings.theta.frequencyBand; 
                                        % Get the theta phase for each spike
                                        thetaData = getPhaseLocking(dataForPhaseLocking, thetaBand);
                                        allTrialPhases{iField}{iTrial} = thetaData.allPhs';

                                    elseif isempty(spikesByTrial) == 1; 
                                        allTrialPhases{iField}{iTrial} = []; 
                                    end
                                    
                                    % Get theta modulation for bursts
                                    if isempty(burstsByTrial) == 0; 
                                        % Define the rest of the inputs for
                                        % the theta modulation script
                                        dataForPhaseLocking.spikeTimes = round(burstsByTrial);
                                        thetaBand = settings.theta.frequencyBand; 
                                        % Get the theta phase for each spike
                                        thetaData = getPhaseLocking(dataForPhaseLocking, thetaBand);
                                        allTrialPhases_bursts{iField}{iTrial} = thetaData.allPhs';

                                    elseif isempty(burstsByTrial) == 1; 
                                        allTrialPhases_bursts{iField}{iTrial} = []; 
                                    end
                                    
                                    % Get theta modulation for singles
                                    if isempty(singlesByTrial) == 0; 
                                        % Define the rest of the inputs for
                                        % the theta modulation script
                                        dataForPhaseLocking.spikeTimes = round(singlesByTrial);
                                        thetaBand = settings.theta.frequencyBand; 
                                        % Get the theta phase for each spike
                                        thetaData = getPhaseLocking(dataForPhaseLocking, thetaBand);
                                        allTrialPhases_singles{iField}{iTrial} = thetaData.allPhs';

                                    elseif isempty(singlesByTrial) == 1; 
                                        allTrialPhases_singles{iField}{iTrial} = []; 
                                    end
                                end

                                % Save the outputs for all spikes, bursts, and singles
                                if strcmp(directions(iDir), 'cw') == 1;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.phases.cw = allTrialPhases;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.burstsPhases.cw = allTrialPhases_bursts;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.singlesPhases.cw = allTrialPhases_singles;
                                elseif strcmp(directions(iDir), 'ccw') == 1; 
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.phases.ccw = allTrialPhases;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.burstsPhases.ccw = allTrialPhases_bursts;
                                    data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.singlesPhases.ccw = allTrialPhases_singles;
                                end
                            
                                % Get the all-trials theta modulation
                                if isempty(allSpikes_allTrials) == 0; 
                                    % Define the rest of the inputs for the theta modulation script
                                    dataForPhaseLocking_allTrials.spikeTimes = round(allSpikes_allTrials);
                                    dataForPhaseLocking_allTrials.LFPsamples = dataForPhaseLocking.LFPsamples;
                                    dataForPhaseLocking_allTrials.LFPtimes = dataForPhaseLocking.LFPtimes;
                                    dataForPhaseLocking_allTrials.samplingFreq = dataForPhaseLocking.samplingFreq; 
                                    thetaBand = settings.theta.frequencyBand; 

                                    % Get the theta locking for all spikes for that field
                                    thetaData_allTrials = getPhaseLocking(dataForPhaseLocking_allTrials, thetaBand);

                                    % Save the data by field
                                    if strcmp(directions(iDir), 'cw') == 1;
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.mvl.cw = thetaData_allTrials.mvl;
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.dir.cw = thetaData_allTrials.dir;
                                    elseif strcmp(directions(iDir), 'ccw') == 1; 
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.mvl.ccw = thetaData_allTrials.mvl;
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.dir.ccw = thetaData_allTrials.dir;
                                    end

                                    % Save the data for the population
                                    populationMVL = [populationMVL, thetaData_allTrials.mvl]; 
                                    populationDir = [populationDir, thetaData_allTrials.dir];

                                    % Get the data for the population rose plot
                                    fieldHistCount = (hist(cell2mat(allTrialPhases{iField}), settings.theta.numBins))/length(cell2mat(allTrialPhases{iField}));
                                    populationHistogram = [populationHistogram; fieldHistCount];
                                end
                                
                                % Get the all-trials bursts theta modulation
                                if isempty(allBursts_allTrials) == 0; 
                                    % Define the rest of the inputs for the theta modulation script
                                    dataForPhaseLocking_allTrials.spikeTimes = round(allBursts_allTrials);
                                    dataForPhaseLocking_allTrials.LFPsamples = dataForPhaseLocking.LFPsamples;
                                    dataForPhaseLocking_allTrials.LFPtimes = dataForPhaseLocking.LFPtimes;
                                    dataForPhaseLocking_allTrials.samplingFreq = dataForPhaseLocking.samplingFreq; 
                                    thetaBand = settings.theta.frequencyBand; 

                                    % Get the theta locking for all spikes for that field
                                    thetaData_allTrials = getPhaseLocking(dataForPhaseLocking_allTrials, thetaBand);

                                    % Save the data by field
                                    if strcmp(directions(iDir), 'cw') == 1;
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.burstsMVL.cw = thetaData_allTrials.mvl;
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.burstsDir.cw = thetaData_allTrials.dir;
                                    elseif strcmp(directions(iDir), 'ccw') == 1; 
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.burstsMVL.ccw = thetaData_allTrials.mvl;
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.burstsDir.ccw = thetaData_allTrials.dir;
                                    end

                                    % Save the data for the population
                                    populationBurstsMVL = [populationBurstsMVL, thetaData_allTrials.mvl]; 
                                    populationBurstsDir = [populationBurstsDir, thetaData_allTrials.dir];

                                    % Get the data for the population rose plot
                                    fieldHistCount = (hist(cell2mat(allTrialPhases_bursts{iField}), settings.theta.numBins))/length(cell2mat(allTrialPhases_bursts{iField}));
                                    populationBurstsHistogram = [populationBurstsHistogram; fieldHistCount];
                                end
                                
                                % Get the all-trials singles theta modulation
                                if isempty(allSingles_allTrials) == 0; 
                                    % Define the rest of the inputs for the theta modulation script
                                    dataForPhaseLocking_allTrials.spikeTimes = round(allSingles_allTrials);
                                    dataForPhaseLocking_allTrials.LFPsamples = dataForPhaseLocking.LFPsamples;
                                    dataForPhaseLocking_allTrials.LFPtimes = dataForPhaseLocking.LFPtimes;
                                    dataForPhaseLocking_allTrials.samplingFreq = dataForPhaseLocking.samplingFreq; 
                                    thetaBand = settings.theta.frequencyBand; 

                                    % Get the theta locking for all spikes for that field
                                    thetaData_allTrials = getPhaseLocking(dataForPhaseLocking_allTrials, thetaBand);

                                    % Save the data by field
                                    if strcmp(directions(iDir), 'cw') == 1;
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.singlesMVL.cw = thetaData_allTrials.mvl;
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.singlesDir.cw = thetaData_allTrials.dir;
                                    elseif strcmp(directions(iDir), 'ccw') == 1; 
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.singlesMVL.ccw = thetaData_allTrials.mvl;
                                        data.cellData.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).theta.singlesDir.ccw = thetaData_allTrials.dir;
                                    end

                                    % Save the data for the population
                                    populationSinglesMVL = [populationSinglesMVL, thetaData_allTrials.mvl]; 
                                    populationSinglesDir = [populationSinglesDir, thetaData_allTrials.dir];

                                    % Get the data for the population rose plot
                                    fieldHistCount = (hist(cell2mat(allTrialPhases_singles{iField}), settings.theta.numBins))/length(cell2mat(allTrialPhases_singles{iField}));
                                    populationSinglesHistogram = [populationSinglesHistogram; fieldHistCount];
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % Get the average MVL represented as a rose plot across the population
        binnedPhases = -pi+pi/12 : pi/12 : pi;
        binIndices = discretize(deg2rad(populationDir), binnedPhases);

        %[counts, binIndices] = histcounts(deg2rad(populationDir), binnedPhases); 
        averageMVLrosePlot = NaN(1, length(binnedPhases)); 
        for iBins = 1:length(binIndices);
            check1 = averageMVLrosePlot(binIndices(iBins))
            check2 = populationMVL(iBins)
            averageMVLrosePlot(binIndices(iBins)) = nanmean([averageMVLrosePlot(binIndices(iBins)), populationMVL(iBins)]);
            display(['finished ', num2str(iBins), ' bin of ', num2str(length(binIndices))])
        end
        
        % Save the population data
        data.populationData(iGenotype).MVL = populationMVL;
        data.populationData(iGenotype).preferredDirection = populationDir;
        data.populationData(iGenotype).popAveRosePlot = averageMVLrosePlot;
        data.populationData(iGenotype).burstsMVL = populationBurstsMVL;
        data.populationData(iGenotype).burstsPreferredDir = populationBurstsDir;
        data.populationData(iGenotype).singlesMVL = populationSinglesMVL;
        data.populationData(iGenotype).singlesPreferredDir = populationSinglesDir;
    end
    
    %% Step 2: Save
    saveFile_v1_20240718(processedDataPath, data, settings, 'theta') 
    
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

    function numField = whichField(fieldsToAnalyze, spikes)  
        % Loop through fields
        if strcmp(fieldsToAnalyze, 'all fields') == 1;
            numField = length(spikes); 
        elseif strcmp(fieldsToAnalyze, 'best field') == 1;
            numField = 1; 
        end
    end

end
    
    
    