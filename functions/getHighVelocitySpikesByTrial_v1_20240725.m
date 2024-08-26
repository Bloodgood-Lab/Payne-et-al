function data = getHighVelocitySpikesByTrial_v1_20240725(data, settings, processedDataPath)
    % Gets the spike times that occur at high velocities and appends that
    % to the data structure
    % Written by Anja Payne
    % Last Modified: 07/30/2024

    % Inputs:
    %   1) data: the matlab structure where the spike times are stored
    %   2) settings: the matlab settings file where the velocity thresholds 
    %      are stored
    
    % Outputs:
    %   1) data: the structure with the high-velocity spike times included

    % Steps:
    %   1) Loop through genotypes, animals, and cells to get the velocity
    %   2) Use the thresholds provided by the user in the settings file to
    %      get the spike times that occur during high-velocity
    %   3) Split the spike times into trials
    %   4) For each trial, determine if it's clockwise or counter-clockwise
    %   5) Ask the user if they want to save the newly generated data
    
    stepSize = settings.velocity.samplingRate * settings.velocity.timeToAverage; % step size is determined by number of frames per second and the time to average over
    directoryPath = '';
    for iGenotype = 1:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype}); 
        for iAnimal = 1:length(genotypeData); 
            tempAnimalStruct = struct(); 
            if isempty(genotypeData{iAnimal}) == 1; 
                continue
            else
                [~,n] = size(genotypeData{iAnimal});
                for iCluster = 1:n; 
                    %% Step 1: Get the velocity
                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                    newDirectoryPath = genotypeData{iAnimal}(iCluster).directory{1};
                    if strcmp(directoryPath, newDirectoryPath) == 1;
                        % The velocity for this session has already been
                        % calculated, do nothing
                        newSession = 0; % not a new session
                    elseif strcmp(directoryPath, newDirectoryPath) == 0;
                        newSession = 1; % it is a new session
                        clear x_pos_cm y_pos_cm t_clipped; 
                        
                        directoryPath = newDirectoryPath;

                        % Load the position data for that session
                        idcs = strfind(directoryPath, '\');
                        parent_directory = directoryPath(1:idcs(end)-1);
                        positionStruct = load([parent_directory, '\Position.mat']);

                        % Exclude any position data that occurs before the video turns on.
                        clipped_indices = find(positionStruct.aligned_pos(:,1) ~= -1000);
                        x_clipped = positionStruct.aligned_pos(clipped_indices, 2);
                        y_clipped = positionStruct.aligned_pos(clipped_indices, 1); 
                        t_clipped = positionStruct.aligned_pos_ts(clipped_indices);
                        clear positionStruct; % clear variable to free up space
                        
                        % Center the track at coordinates 0,0 and convert
                        % to cm
                        [x_pos_centered, y_pos_centered] = centerBox(x_clipped, y_clipped); 
                        clear x_clipped y_clipped; % clear variables to free up space
                        widthConversionFactor = settings.rateMaps.trackWidth/(max(x_pos_centered) - min(x_pos_centered)); 
                        lengthConversionFactor = settings.rateMaps.trackLength/(max(y_pos_centered) - min(y_pos_centered)); 
                        x_pos_cm = x_pos_centered*widthConversionFactor;
                        y_pos_cm = y_pos_centered*lengthConversionFactor;
                        clear x_pos_centered y_pos_centered; % clear variables to free up space; 
                        
                        % Get the velocity and identify high-velocity times
                        velocity = getVelocity(x_pos_cm, y_pos_cm, stepSize); 
                        indxHighVelocity = velocity > settings.velocity.threshold; 
                        timeHighVelocity = t_clipped(indxHighVelocity == 1);
                    end

                    %% Step 2: Get spikes that occur during running
                    runningData = struct(); 
                    cellSpikeTimes = genotypeData{iAnimal}(iCluster).spikeTimes;
                    inter_data = intersect(round(cellSpikeTimes), round(timeHighVelocity)); 
                    inter_ind = find(ismember(round(cellSpikeTimes), inter_data)); 
                    runningData.highVelocitySpikes = cellSpikeTimes(inter_ind);
                    runningData.allSpikes = cellSpikeTimes;
                    runningData.x = x_pos_cm; 
                    runningData.y = y_pos_cm; 
                    runningData.t = t_clipped; 
                    
                    %% Step 3: Split the spikes into trials
                    splitByTrialsData = splitFilesByTrials(runningData); 
                    % Organize one structure to include the high velocity
                    % data and one structure to include all data
                    splitByTrials_highVelocity.s = splitByTrialsData.highVelSpks;
                    splitByTrials_highVelocity.x = splitByTrialsData.x;
                    splitByTrials_highVelocity.y = splitByTrialsData.y;
                    splitByTrials_highVelocity.t = splitByTrialsData.t;
                    splitByTrials_allData.s = splitByTrialsData.allSpks;
                    splitByTrials_allData.x = splitByTrialsData.x;
                    splitByTrials_allData.y = splitByTrialsData.y;
                    splitByTrials_allData.t = splitByTrialsData.t;  
                    clear splitByTrialsData;

                    %% Step 4: Split into clockwise and counterclockwise
                    bins = settings.rateMaps.binSize; 
                    binnedSpikesByDirection_highVelocity = splitCWandCCWandBinSpikes(splitByTrials_highVelocity, bins);
                    binnedSpikesByDirection_allSpikes = splitCWandCCWandBinSpikes(splitByTrials_allData, bins);
                    clear splitByTrials_highVelocity splitByTrials_allData
                    
                    % Append the trials
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).trials.cw = binnedSpikesByDirection_allSpikes.trialNumber.cw;
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).trials.ccw = binnedSpikesByDirection_allSpikes.trialNumber.ccw; 
                    
                    % Append the high velocity spike times and positions
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).highVelocityData.unbinnedSpkPos.cw = binnedSpikesByDirection_highVelocity.spikeUnbinnedPos.cw;
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).highVelocityData.unbinnedSpkPos.ccw = binnedSpikesByDirection_highVelocity.spikeUnbinnedPos.ccw;
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).highVelocityData.spikePosBins.cw = binnedSpikesByDirection_highVelocity.spikePosBins.cw;
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).highVelocityData.spikePosBins.ccw = binnedSpikesByDirection_highVelocity.spikePosBins.ccw; 
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).highVelocityData.spikeTimesByTrial.cw = binnedSpikesByDirection_highVelocity.spikeTimes.cw; 
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).highVelocityData.spikeTimesByTrial.ccw = binnedSpikesByDirection_highVelocity.spikeTimes.ccw; 
                    
                    % Append spike times and positions at all velocities
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).allVelocities.unbinnedSpkPos.cw = binnedSpikesByDirection_allSpikes.spikeUnbinnedPos.cw;
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).allVelocities.unbinnedSpkPos.ccw = binnedSpikesByDirection_allSpikes.spikeUnbinnedPos.ccw;
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).allVelocities.spikePosBins.cw = binnedSpikesByDirection_allSpikes.spikePosBins.cw;
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).allVelocities.spikePosBins.ccw = binnedSpikesByDirection_allSpikes.spikePosBins.ccw; 
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).allVelocities.spikeTimesByTrial.cw = binnedSpikesByDirection_allSpikes.spikeTimes.cw; 
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).allVelocities.spikeTimesByTrial.ccw = binnedSpikesByDirection_allSpikes.spikeTimes.ccw; 

                    % Since the binned positions are so large, save these
                    % into a separate folder and append the filepath 
                    % Only do this once per session
                    if newSession == 0;
                        % Positions for this session were already saved
                        data.(genotypes{iGenotype}){iAnimal}(iCluster).posBinFile = posSaveFile; 
                    elseif newSession == 1;
                        saveFolder = [processedDataPath, '\binnedPositions'];
                        if ~exist(saveFolder, 'dir')
                            % Create the folder if it does not exist
                            mkdir(saveFolder);
                        end
                        binnedPosition.cw = binnedSpikesByDirection_allSpikes.posBins.cw;
                        binnedPosition.ccw = binnedSpikesByDirection_allSpikes.posBins.ccw;
                        clear binnedSpikesByDirection_highVelocity binnedSpikesByDirection_allSpikes
                        
                        tempDirectoryPath = fileparts(newDirectoryPath); 
                        [~, saveFolderName, ~] = fileparts(tempDirectoryPath);
                        posSaveFile = [saveFolder, '\', saveFolderName, '_binnedPos'];
                        save(posSaveFile, 'binnedPosition'); 
                        data.(genotypes{iGenotype}){iAnimal}(iCluster).posBinFile = posSaveFile; 
                    end
                end
            end
        end
    end
    
    %% Step 2: Step 2: Save
    saveFile_v1_20240718(processedDataPath, data, settings, 'highVelocitySpikeTimes') 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [posx_center, posy_center] = centerBox(posx, posy)
        % Script to find the center of the behavioral box
        % Written in 2013 by members of Stefan Leutgeb's lab
        % Inputs:
        %   1) posx: x-position in pixel units
        %   2) posy: y-position in pixel units
        % Outputs:
        %   1) posx_center: centered position in cm
        %   2) posy_center: centered position in cm

        % Find border values for box for session 1
        maxX = max(posx);
        minX = min(posx);
        maxY = max(posy);
        minY = min(posy);

        % Set the corners of the reference box
        NE = [maxX, maxY];
        NW = [minX, maxY];
        SW = [minX, minY];
        SE = [maxX, minY];

        % The centre will be at the point of interception by the corner diagonals
        a = (NE(2)-SW(2))/(NE(1)-SW(1)); % Slope for the NE-SW diagonal
        b = (SE(2)-NW(2))/(SE(1)-NW(1)); % Slope for the SE-NW diagonal
        c = SW(2);
        d = NW(2);
        x = (d-c+a*SW(1)-b*NW(1))/(a-b); % X-coord of centre
        y = a*(x-SW(1))+c; % Y-coord of centre
        center = [x,y];

        % Center both boxes according to the coordinates to the first box
        posx_center = posx - center(1);
        posy_center = posy - center(2);
    end

    function velocity = getVelocity(x, y, stepSize)
        % Given the x and y position and the stepSize, get the velocity
        % Inputs:
        %   1) x: x-position in cm
        %   2) y: y-position in cm
        % Outputs:
        %   1) velocity: velocity in cm/sec sampled to match the position
        %      data

        % Calculate the velocity for that session
        x = x + min(x); % Add the min to ensure only positive values
        y = y + min(y); % Add the min to ensure only positive values
        velocity = NaN(1, length(x)); 
        % Step through using the defined stepSize and determine the velocity for
        % each frame
        for iWindow = (1+stepSize/2):(length(x)-(stepSize/2)-1);
            x_diff = abs(x(iWindow - (stepSize/2)) - x(iWindow + (stepSize/2)));
            y_diff = abs(y(iWindow - (stepSize/2)) - y(iWindow + (stepSize/2)));
            velocity(iWindow) = (sqrt(x_diff.^2 + y_diff.^2));
        end
        
    end
    
    function output = splitFilesByTrials(input)
        % Given spike times and position, split the spike times into trials
        % Inputs:
        %   1) allSpikes: spike times
        %   2) allX: x-positions
        %   3) allY: y-positions
        %   4) allT: all position timepoints
        % Outputs:
        %   1) Spikes as a matrix with dimensions number of trials x time
        %   2) X-position as a matrix with dimensions number of trials x time
        %   3) Y-position as a matrix with dimensions number of trials x time
        %   4) Time as a matrix with dimensions number of trials x time
        % Steps:
        %   1) Find putative 'trial starts' by finding every time the animal's 
        %      x-position is between 0 and |10| and the animal's y-position is between 
        %      -30 and -40. 
        %   2) Loop through the putative trial starts and, beginning with the first
        %      one, narrow the trial start selections by ignoring the next 20
        %      seconds. 
        %   3) Store that trial in one row of a matrix. 
        %   4) For each trial, make sure that the animal also passes through the
        %      point that is at least greater than 30. If it doesn't, exclude this
        %      trial. 
        %   5) Save the animal's x-position data, y-position data, and time as a 
        %      matrix with the y-axis being trial number. 

        output = struct(); 
        highVspks = input.highVelocitySpikes;
        allSpikes = input.allSpikes; 
        allX = input.x;
        allY = input.y; 
        allT = input.t;
        
        % Step 1: Find the putative trial start times
        xValues = find(abs(allX)<10); 
        yValues1 = find(allY > -40);
        yValues2 = find(allY < -30);
        yValues = intersect(yValues1, yValues2);
        putativeTrialStart = intersect(xValues, yValues); 

        % Step 2: Refine the trial start times
        trialJump = 50000; % ~3 seconds; used historically
        count = 1; 
        trialStart = NaN(1, 10000); 
        
        for itrialStart = 1:length(putativeTrialStart)-1; 
            % Find the points where the animal's position 'jumps.' That is, when
            % a new trial starts. 
            tempJump = diff([putativeTrialStart(itrialStart), putativeTrialStart(itrialStart+1)]);
            if tempJump > trialJump;
                trialStart(count) = putativeTrialStart(itrialStart); 
                count = count +1; 
            end
        end

        % Step 3: Store that trial in one row of a matrix
        % Save each trial in its own row of a matrix
        trialStart(isnan(trialStart)) = []; 
        count2 = 1; allTrialsY = {}; allTrialsX = {}; 
        allTrialsSpikes = {}; allTrialsT = {};
        for iTrials = 1:length(trialStart)-1; 
            tempTrialY = allY(trialStart(iTrials):trialStart(iTrials+1)-1); 
            % Step 4: Make sure that the trial includes a timepoint where the
            % animal runs across the back of the track
            if sum(tempTrialY > 30) ~= 0;
                allTrialsY{count2} = tempTrialY; 
                tempTrialX = allX(trialStart(iTrials):trialStart(iTrials+1)-1);
                allTrialsX{count2} = tempTrialX; 
                tempTrialT = allT(trialStart(iTrials):trialStart(iTrials+1)-1); 
                allTrialsT{count2} = tempTrialT; 
                tempTrialSpikes = intersect(round(allSpikes), round(tempTrialT)); 
                allTrialsSpikes{count2} = tempTrialSpikes;
                tempHighVelSpikes = intersect(round(highVspks), round(tempTrialT)); 
                highVelSpikes{count2} = tempHighVelSpikes;
                count2 = count2 + 1;
            end
        end

        output.allSpks = allTrialsSpikes;
        output.highVelSpks = highVelSpikes;
        output.x = allTrialsX; 
        output.y = allTrialsY;
        output.t = allTrialsT; 

    end
    
    function output = splitCWandCCWandBinSpikes(input, bins)
        % Splits trials in to clockwise and counterclockwise
        % Inputs: 
        %   1) input.spikes = spike times
        %   2) input.x = x-position
        %   3) input.y = y-position
        %   4) input.t = position timepoints
        % Outputs:
        %   1) output.spikeUnbinnedPos: the unbinned position for every
        %      spike
        %   2) output.posBins: the binned position value for every
        %      position, dimensions will match the position dimensions from 
        %      the video 
        %   3) output.spikePosBins: the binned position for every spike,
        %      dimensions will match the number of spikes
        %   4) output.trialNumber: the corresponding trial numbers
        
        sTrials = input.s;
        xTrials = input.x;
        yTrials = input.y;
        tTrials = input.t;
                
        count_cw = 1; count_ccw = 1; 
        cw_posBins = {}; ccw_posBins = {}; cw_spikePosBins = {}; ccw_spikePosBins = {}; 
        cw_spikes = []; ccw_spikes = []; cw_trials = []; ccw_trials = []; 
        for iTrial = 1:length(xTrials); 
            
            % Define all points that are not equal to zero or to NaN
            pts = isfinite(xTrials{iTrial}(:)) & isfinite(yTrials{iTrial}(:)) & ~(xTrials{iTrial}(:)==0 & yTrials{iTrial}(:)==0);
            
            % Calculate the angle such that the bottom of the map (where the reward
            % is) is -pi and the position increases as the animal moves
            % counter-clockwise. 
            % If you were plotting x by y to show the animal's two dimensional
            % trajectory then plotting y by -x would give you a map that is rotated
            % -90 degrees such that the start point is now on the left instead of
            % at the bottom. This will map onto the unit circle such that the map 
            % will go from -pi to pi. So for atan2, y should be -x{iTrials} and x 
            % should be y{iTrials}. 
            pts_angle = atan2(-xTrials{iTrial}(pts), yTrials{iTrial}(pts)); 

            % Bin the angular positions into 4cm bins
            numberBins = 264/bins; 
            [~, ~, posBins{iTrial}] = histcounts(pts_angle, linspace(-pi, pi, numberBins)); 
            
            % Find the x and y position of every spike 
            [~, spike_ind] = intersect(round(tTrials{iTrial}), sTrials{iTrial});
            xSpike = xTrials{iTrial}(spike_ind); 
            ySpike = yTrials{iTrial}(spike_ind);
            % Create an array that is the length of the time array where each spike
            % point is denoted as 1 and the rest is 0. 
            spike_array = zeros(1, length(tTrials{iTrial}));
            spike_array(spike_ind) = 1; 
            % Only keep the values that are finite (not NaN or 0)
            spike_finitePts = isfinite(xSpike(:)) & isfinite(ySpike(:)) & ~(xSpike(:)==0 & ySpike(:)==0);
            spike_posLogical = spike_array(pts);
            % For each spike, find the position time array index
            spike_pos_pts{iTrial} = find(spike_posLogical == 1);
            spike_pts{iTrial} = find(spike_finitePts == 1);
            include_spikes = sTrials{iTrial}(spike_pts{iTrial});
            
            % Find the bin associated with each point (10 msec) of time
            pts_pos = isfinite(xTrials{iTrial}) & isfinite(yTrials{iTrial}) & ~(xTrials{iTrial}==0 & yTrials{iTrial} == 0); 
            pts_pos_angle = atan2(-xTrials{iTrial}(pts_pos), yTrials{iTrial}(pts_pos));
            [~,~,pos_bins_temp] = histcounts(pts_pos_angle, linspace(-pi, pi, numberBins)); 
            
            % Split position_bin into the two directions
            % Check the position on the middle part of the track to see if
            % the animal is running clockwise or counter-clockwise
            testArray = posBins{iTrial}(posBins{iTrial} > (0.25*numberBins) & posBins{iTrial} < (0.75*numberBins));
            tempSum = testArray(end) - testArray(1);       
            if tempSum < 0; % if the values are decreasing the animal is running clockwise
               cw_unbinnedSpikePos{count_cw} = rad2deg(pts_angle(spike_pos_pts{iTrial})); 
               cw_posBins{count_cw} = posBins{iTrial}; 
               cw_spikePosBins{count_cw} = posBins{iTrial}(spike_pos_pts{iTrial}); 
               cw_spikeTimes{count_cw} = include_spikes;
               cw_trials(count_cw) = iTrial; 
               count_cw = count_cw + 1;
            elseif tempSum >= 0; % if the values are increasing the animal is running counter-clockwise
               ccw_unbinnedSpikePos{count_ccw} = rad2deg(pts_angle(spike_pos_pts{iTrial}));
               ccw_posBins{count_ccw} = posBins{iTrial}; 
               ccw_spikePosBins{count_ccw} = posBins{iTrial}(spike_pos_pts{iTrial}); 
               ccw_spikeTimes{count_ccw} = include_spikes;
               ccw_trials(count_ccw) = iTrial; 
               count_ccw = count_ccw + 1;
            end
        end
        
        output.spikeUnbinnedPos.cw = cw_unbinnedSpikePos; 
        output.spikeUnbinnedPos.ccw = ccw_unbinnedSpikePos; 
        output.posBins.cw = cw_posBins;
        output.posBins.ccw = ccw_posBins; 
        output.spikePosBins.cw = cw_spikePosBins; 
        output.spikePosBins.ccw = ccw_spikePosBins; 
        output.spikeTimes.cw = cw_spikeTimes;
        output.spikeTimes.ccw = ccw_spikeTimes;
        output.trialNumber.cw = cw_trials; 
        output.trialNumber.ccw = ccw_trials; 

    end
end
