function data = getInFieldSpikes_v1_20240805(data, settings)
    % Gets the in-field spikes
    % Written by Anja Payne
    % Last Modified: 08/05/2024

    % Inputs:
    %   1) data: the matlab structure where the spikes-by-trial and the
    %      barcode are saved
    
    % Outputs:
    %   1) data: the structure with in-field spikes appended

    % Steps:
    %   1) Loop through genotypes, animals, and cells and extract the 
    %      in-field spikes
    %   2) Ask the user if they want to save the newly generated data
    
    %% Step 1: Get the in-field spikes
    for iGenotype = 1:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype}); 
        
        % Only run for high-firing cells
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
                        continue
                    else
                        display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                        directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode);
                        for iDir = 1:length(directions);
                            if strcmp(directions(iDir), 'cw') == 1; 
                                barcode = FRdata{iAnimal}(iCluster).spatialMetrics.barcode.cw;
                                spkPos = FRdata{iAnimal}(iCluster).binnedSpikesByTrial.allVelocities.binnedSpkPos.cw;
                                spkTimes = FRdata{iAnimal}(iCluster).binnedSpikesByTrial.allVelocities.binnedSpkTimes.cw;
                            elseif strcmp(directions(iDir), 'ccw') == 1; 
                                barcode = FRdata{iAnimal}(iCluster).spatialMetrics.barcode.ccw;
                                spkPos = FRdata{iAnimal}(iCluster).binnedSpikesByTrial.allVelocities.binnedSpkPos.ccw;
                                spkTimes = FRdata{iAnimal}(iCluster).binnedSpikesByTrial.allVelocities.binnedSpkTimes.ccw;
                            end

                            for iField = 1:max(barcode);
                                inFieldSpkTimes = {}; 
                                inFieldIndex = find(barcode == iField); 
                                for iTrial = 1:length(spkPos); 
                                    inFieldLogical = ismember(spkPos{iTrial}, inFieldIndex); 
                                    inFieldSpkTimes{iTrial} = spkTimes{iTrial}(inFieldLogical); 
                                end

                                if strcmp(directions(iDir), 'cw') == 1;
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldSpkTimes.cw{iField} = inFieldSpkTimes;
                                elseif strcmp(directions(iDir), 'ccw') == 1; 
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldSpkTimes.ccw{iField} = inFieldSpkTimes;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% Step 2: Save
    saveFile_v1_20240718(data, settings, 'inFieldSpkTimes') 
    