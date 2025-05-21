function data = getSpikeTimingData_v1_20241124(data, settings, processedDataPath)
    % Gets the spike timing data (bursts vs singles)
    % Written by Anja Payne
    % Last Modified: 11/22/2024

    % Inputs:
    %   1) data: the matlab structure where the spike times are saved
    
    % Outputs:
    %   1) data: the structure with bursting and single data saved

    % Steps:
    %   1) Loop through genotypes, animals, and cells and get the mean and
    %      max firing rates over time
    %   2) Extract the bursts and singles
    %   3) Ask the user if they want to save the newly generated data
    
    %%% Step 1: Get the mean and max firing rates over time
    for iGenotype = 1:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype}); 
        for iAnimal = 1:length(genotypeData); 
            % Skip if empty
            if isempty(genotypeData{iAnimal}) == 1; 
                continue
            else
                [~,n] = size(genotypeData{iAnimal});
                for iCluster = 1:n;
                    display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                    directions = {'cw', 'ccw'}; 
                    for iDir = 1:length(directions);
                        outputData = assignVariableByDirection_v1_20240905(genotypeData{iAnimal}(iCluster), directions(iDir), 'spikeTiming');
                        spkTimes = outputData.spkTimes;

                        
                        
                        
                        singles = cell(1, length(spkTimes)); bursts = cell(1, length(spkTimes)); 
                        for iTrial = 1:length(spkTimes); 
                            if isempty(spkTimes{iTrial}) == 0 %if there are spikes
                                
                                %%% Step 2: Extract bursts and singles
                                singles{iTrial} = NaN(1, length(spkTimes{iTrial}));  
                                bursts{iTrial} = NaN(1, length(spkTimes{iTrial})); 
                                if length(spkTimes{iTrial}) == 1 %if there's only one element it's a single
                                    singles{iTrial}(1) = spkTimes{iTrial}(1);
                                else
                                    spike_ISIs = diff(spkTimes{iTrial}); %get the ISIs
                                    spike_ISIs = [NaN; spike_ISIs]; %pad it so that it's the same length as spkTimes
                                    for spikes = 2:length(spike_ISIs);
                                        % If the ISI is greater than 10 then count
                                        % it as a single. If the ISI is less than 
                                        % 10 then count it as a burst
                                        if spike_ISIs(spikes) > 10;% && spike_ISIs(spikes+1) > 10;
                                           singles{iTrial}(spikes) = spkTimes{iTrial}(spikes);  
                                           if spikes == 2; 
                                               singles{iTrial}(1) = spkTimes{iTrial}(1); 
                                           end
                                        end
                                        if spike_ISIs(spikes) <= 10;
                                           bursts{iTrial}(spikes) = spkTimes{iTrial}(spikes); 
                                           bursts{iTrial}(spikes-1) = spkTimes{iTrial}(spikes-1); 
                                           singles{iTrial}(spikes-1) = NaN; 
                                        end
                                    end
                                end
                            end
                        end
                        if strcmp(directions(iDir), 'cw') == 1;
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).spikeTiming.bursts.cw = bursts;
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).spikeTiming.singles.cw = singles;
                        elseif strcmp(directions(iDir), 'ccw') == 1; 
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).spikeTiming.bursts.ccw = bursts;
                            data.(genotypes{iGenotype}){iAnimal}(iCluster).spikeTiming.singles.ccw = singles;
                        end
                    end
                end
            end
        end
    end
    
    %% Step 2: Save
    saveFile_v1_20240718(processedDataPath, data, settings, 'spikeTiming'); 
    