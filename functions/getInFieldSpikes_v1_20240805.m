function data = getInFieldSpikes_v1_20240805(data, settings, processedDataPath)
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
                        directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original);
                        for iDir = 1:length(directions);
                            outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir), 'getInField');
                            barcode = outputData.barcode; binnedSpkPos = outputData.binnedSpkPosForInField; 
                            spkPos = outputData.spkPosForInField; spkTimes = outputData.spkTimes;
                            bursts = outputData.bursts; singles = outputData.singles; 
                            
                            for iField = 1:max(barcode);
                                inFieldLogical = []; inFieldSpkPos = {}; inFieldSpkTimes = {}; inFieldSpkPosBinned = {}; 
                                inFieldIndex = find(barcode == iField); 
                                for iTrial = 1:length(binnedSpkPos); 
                                    inFieldLogical = ismember(binnedSpkPos{iTrial}, inFieldIndex); 
                                    inFieldSpkTimes{iTrial} = spkTimes{iTrial}(inFieldLogical);
                                    inFieldSpkPosBinned{iTrial} = binnedSpkPos{iTrial}(inFieldLogical); 
                                    inFieldSpkPos{iTrial} = spkPos{iTrial}(inFieldLogical); 
                                    inFieldBursts{iTrial} = bursts{iTrial}(inFieldLogical); 
                                    inFieldSingles{iTrial} = singles{iTrial}(inFieldLogical); 
                                end

                                if strcmp(directions(iDir), 'cw') == 1;
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldSpkTimes.cw{iField} = inFieldSpkTimes;
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldBinnedSpkPos.cw{iField} = inFieldSpkPosBinned;
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldSpkPos.cw{iField} = inFieldSpkPos;
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldBursts.cw{iField} = inFieldBursts;
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldSingles.cw{iField} = inFieldSingles;
                                elseif strcmp(directions(iDir), 'ccw') == 1; 
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldSpkTimes.ccw{iField} = inFieldSpkTimes;
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldBinnedSpkPos.ccw{iField} = inFieldSpkPosBinned;
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldSpkPos.ccw{iField} = inFieldSpkPos;
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldBursts.ccw{iField} = inFieldBursts;
                                    data.(genotypes{iGenotype}).highFiring{iAnimal}(iCluster).inField.inFieldSingles.ccw{iField} = inFieldSingles;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% Step 2: Save
    saveFile_v1_20240718(processedDataPath, data, settings, 'inFieldSpkTimes') 
    