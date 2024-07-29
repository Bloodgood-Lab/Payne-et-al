function data = getHighVelocitySpikesByTrial_v1_20240725(data, settings)
    % Gets the spike times that occur at high velocities and appends that
    % to the data structure
    % Written by Anja Payne
    % Last Modified: 07/25/2024

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
    
    stepSize = 16000 * settings.velocity.timeToAverage; % step size is determined by number of frames per second and the time to average over
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
                    %% Step 1: Get the velocity
                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);

                    newDirectoryPath = genotypeData{iAnimal}(iCluster).directory{1};
                    if strcmp(directoryPath, newDirectoryPath) == 1;
                        % The velocity for this session has already been
                        % calculated, do nothing
                    elseif strcmp(directoryPath, newDirectoryPath) == 0;
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
                        
                        % Center the track at coordinates 0,0 and convert
                        % to cm
                        [x_pos_centered, y_pos_centered] = centerBox(x_clipped, y_clipped);
                        widthConversionFactor = settings.rateMaps.trackWidth/(max(x_pos_centered) - min(x_pos_centered)); 
                        lengthConversionFactor = settings.rateMaps.trackLength/(max(y_pos_centered) - min(y_pos_centered)); 
                        x_pos_cm = x_pos_centered*widthConversionFactor;
                        y_pos_cm = y_pos_centered*lengthConversionFactor;
                        
                        % Get the velocity and identiy high-velocity times
                        velocity = getVelocity(x_pos_cm, y_pos_cm, stepSize); 
                        indxHighVelocity = velocity > settings.velocity.threshold; 
                        timeHighVelocity = t_clipped(indxHighVelocity == 1);
                    end

                    %% Step 2: Get spikes that occur during running
                    highVelocityData = struct(); 
                    cellSpikeTimes = genotypeData{iAnimal}(iCluster).spikeTimes;
                    inter_data = intersect(round(cellSpikeTimes), round(timeHighVelocity)); 
                    inter_ind = find(ismember(round(cellSpikeTimes), inter_data)); 
                    highVelocityData.spikes = cellSpikeTimes(inter_ind);
                    highVelocityData.x = x_pos_cm;
                    highVelocityData.y = y_pos_cm;
                    highVelocityData.t = t_clipped; 
                    

                    %% Step 3: Split the spikes into trials
                    splitByTrialsData = splitFilesByTrials(highVelocityData); 
                    
                    %% Step 4: Split into clockwise and counterclockwise
                    bins = settings.rateMaps.binSize; 
                    binnedSpikesByDirection = splitCWandCCWandBinSpikes(splitByTrialsData, bins);

                    % Append to data structure
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).binnedPos.cw = binnedSpikesByDirection.pos.cw;
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).binnedPos.ccw = binnedSpikesByDirection.pos.ccw; 
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).binnedSpikes.cw = binnedSpikesByDirection.spikes_binned.cw;
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).binnedSpikes.ccw = binnedSpikesByDirection.spikes_binned.ccw; 
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).spikesByTrial.cw = binnedSpikesByDirection.spikes.cw; 
                    data.(genotypes{iGenotype}){iAnimal}(iCluster).spikesByTrial.ccw = binnedSpikesByDirection.spikes.ccw;
                    
                end
            end
        end
    end
    
    
    
    %% Step 2: Step 2: Save
    saveFile_v1_20240718(data, settings, 'highVelocitySpikeTimes') 
    
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
        allSpikes = input.spikes;
        allX = input.x;
        allY = input.y; 
        allT = input.t;
        
        % Step 1: Find the putative trial start times
        xValues = find(abs(allX)<10); 
        yValues1 = find(allY > -40);
        yValues2 = find(allY < -30);
        yValues = intersect(yValues1, yValues2);
        putativeTrialStart = intersect(xValues, yValues); % TRACED ISSUE BACK TO HERE (MAYBE EARLIER)

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
                count2 = count2 + 1;
            end
        end

        output.spikes = allTrialsSpikes;
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
        %   1) Find the bin location of each spike on the linearized track
        %   2) Find the bin location of each time point
        %   3) Allocate the data into cw and ccw
        
        sTrials = input.spikes;
        xTrials = input.x;
        yTrials = input.y;
        tTrials = input.t;
        
        count1 = 1; count2 = 1; pos_bins = {}; 
        cw_pos = []; ccw_pos = []; cw_spike_bins = []; ccw_spike_bins = []; cw_spikes = []; ccw_spikes = []; 
        for iTrial = 1:length(xTrials); 
            % Step 1: Find the bin location of each spike on the linearized track
            % Find the x and y position of every spike 
            [~, tSpike] = intersect(round(tTrials{iTrial}), sTrials{iTrial});
            xSpike = xTrials{iTrial}(tSpike);
            ySpike = yTrials{iTrial}(tSpike); 
            % Find the positions that are non NaN            
            pts = isfinite(xSpike(:)) & isfinite(ySpike(:)) & ~(xSpike(:)==0 & ySpike(:)==0); 
            trialSpikes = sTrials{iTrial}(pts); 
            pts_angle = atan2(-xSpike(pts), ySpike(pts));
            numberBins = 264/bins; 
            [~,~,spike_bins] = histcounts(pts_angle, linspace(-pi, pi, numberBins)); 
            
            % Step 2: Find the bin location of each time point
            % Find the bin associated with each point (10 msec) of time
            pts_pos = isfinite(xTrials{iTrial}) & isfinite(yTrials{iTrial}) & ~(xTrials{iTrial}==0 & yTrials{iTrial} == 0); 
            pts_pos_angle = atan2(-xTrials{iTrial}(pts_pos), yTrials{iTrial}(pts_pos));
            [~,~,pos_bins_temp] = histcounts(pts_pos_angle, linspace(-pi, pi, numberBins)); 
            % Subsample the position information so that each bin is 10 msec
            pos_bins{iTrial} = pos_bins_temp(1:160:end); 
            
            % Step 3: Save the data into the two directions
            % Split position_bin into the two directions
            % Check the position on the middle part of the track to see if
            % the animal is running clockwise or counter-clockwise
            testArray = pos_bins{iTrial}(pos_bins{iTrial} > (0.25*numberBins) & pos_bins{iTrial} < (0.75*numberBins));
            tempSum = testArray(end) - testArray(1);                                                               
            if tempSum < 0; % if the values are decreasing the animal is running clockwise
               cw_pos = [cw_pos; pos_bins{iTrial}]; 
               cw_spike_bins = [cw_spike_bins; spike_bins]; 
               cw_spikes = [cw_spikes; trialSpikes];
               count1 = count1 + 1;
            elseif tempSum > 0; % if the values are increasing the animal is running counter-clockwise
               ccw_pos = [ccw_pos; pos_bins{iTrial}]; 
               ccw_spike_bins = [ccw_spike_bins; spike_bins]; 
               ccw_spikes = [ccw_spikes; trialSpikes];
               count2 = count2 + 1;
            end
        end
        
        output.pos.cw = cw_pos; 
        output.pos.ccw = ccw_pos; 
        output.spikes_binned.cw = cw_spike_bins; 
        output.spikes_binned.ccw = ccw_spike_bins; 
        output.spikes.cw = cw_spikes; 
        output.spikes.ccw = ccw_spikes;
        
    end
end
