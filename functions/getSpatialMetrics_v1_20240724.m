function data = getSpatialMetrics_v1_20240724(data, settings, processedDataPath)
    % Gets the spatial metrics from the trial-averaged rate maps
    % Written by Anja Payne
    % Last Modified: 07/30/2024

    % Inputs:
    %   1) data: the matlab structure where the trial-averaged rate
    %      maps are stored
    %   2) settings: the matlab settings file where the thresholds are
    %      stored
    
    % Outputs:
    %   1) data: the structure with spatial metrics (for now just
    %      barcode, place field size, and place field number) appended

    % Steps:
    %   1) Loop through genotypes, animals, and cells and get the barcode,
    %      average size, number of place fields, mean firing rate, and max
    %      firing rate. 
    %   2) Get the max and mean firing rate for each cell
    %   3) Get the spatial information, sparsity, and selectivity
    %   4) Ask the user if they want to save the newly generated data
    
    
    %% Step 1: Get the barcode, average size, and number of place fields
    thresholds(1) = settings.rateMaps.lowThresh; 
    thresholds(2) = settings.rateMaps.highThresh;
    biggerFieldModification = settings.rateMaps.biggerFieldModification;
    smallerFieldModification = settings.rateMaps.smallerFieldModification;
    for iGenotype = 1:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype});
        for iFR = 1:length(fieldnames(genotypeData)); 
            FRoptions = fieldnames(genotypeData); 
            FRdata = genotypeData.(FRoptions{iFR}); 
            for iAnimal = 1:length(FRdata); 
                if isempty(FRdata{iAnimal}) == 1; 
                    continue
                else
                    [~,n] = size(FRdata{iAnimal});
                    for iCluster = 1:n;
                        if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                            continue
                        else
                            display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                            directions = fieldnames(FRdata{iAnimal}(iCluster).rateMaps.rateMap);
                            for iDir = 1:length(directions);
                                % Extract data based on running direction
                                outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir), 'spatialMetrics');
                                map = outputData.map; timeMap = outputData.timeMap; 
                                trialAveragedMap = nanmean(map, 1); 
                                trialAveragedTimeMap = nanmean(timeMap, 1); 
                                                                
                                [barcode, PFsize, PFnumber] = getPlaceFields_v1_20250425(map, thresholds);
                                % Get the position PDF
                                posPDF = trialAveragedTimeMap/nansum(nansum(trialAveragedTimeMap)); 
                                % Get the information and sparsity
                                [info, spars, ~] = getMapStats(trialAveragedMap, posPDF);
                                
                                % For use in later control analyses, modify the 
                                % barcode to make the fields larger or smaller
                                [biggerFieldBarcode, biggerFieldIndices] = modifyPlaceFields(barcode, biggerFieldModification);
                                [smallerFieldBarcode, smallerFieldIndices] = modifyPlaceFields(barcode, smallerFieldModification);
                                
                                % Assign data based on running direction
                                if strcmp(directions(iDir), 'cw') == 1; 
                                    % Define data within the structure
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFsize.cw = PFsize;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFnumber.cw = PFnumber;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.info.cw = info;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.sparsity.cw = spars;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.original.cw = barcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.barcode.cw = biggerFieldBarcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.indices.cw = biggerFieldIndices;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.smaller.barcode.cw = smallerFieldBarcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.smaller.indices.cw = smallerFieldIndices;
                                elseif strcmp(directions(iDir), 'ccw') == 1; 
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFsize.ccw = PFsize;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFnumber.ccw = PFnumber;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.info.ccw = info;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.sparsity.ccw = spars;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.original.ccw = barcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.barcode.ccw = biggerFieldBarcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.indices.ccw = biggerFieldIndices;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.smaller.barcode.ccw = smallerFieldBarcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.smaller.indices.ccw = smallerFieldIndices;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% Step 2: Save
    saveFile_v1_20240718(processedDataPath, data, settings, 'spatialMetrics') 

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [outputBarcode, barcodeIndices] = modifyPlaceFields(inputBarcode, modification)
        % Make the fields either larger or smaller
        
        barcodeIndices = 1:length(inputBarcode); 
        binsToSetToPF = []; PFvalue = [];
        % If the first bin is not zero, set the last bin to that value. Else, if
        % the current bin is zero but the next one isn't, set the current bin to
        % the non-zero value. 
        if modification > 0; % make the field bigger
            for iNumBins = 1:abs(modification);
                startValue = 1; 
                if inputBarcode(1) > 0 && inputBarcode(end) == 0; 
                    binsToSetToPF = [binsToSetToPF, length(inputBarcode)]; 
                    PFvalue = [PFvalue, inputBarcode(1)]; 
                    startValue = 2; 
                end
                for i = startValue:length(inputBarcode)-1;
                    if inputBarcode(i) == 0 && inputBarcode(i+1) > 0; 
                        binsToSetToPF = [binsToSetToPF, i];
                        PFvalue = [PFvalue, inputBarcode(i+1)];
                    end
                end
                inputBarcode(binsToSetToPF) = PFvalue; 
            end

            % If the last bin is not zero, set the first bin to that value. Else, if
            % the current bin is zero but the last one wasn't, set the current bin to
            % the non-zero value. 
            endValue = length(inputBarcode);
            if inputBarcode(end) > 0 && inputBarcode(1) == 0; 
                binsToSetToPF = [binsToSetToPF, 1]; 
                PFvalue = [PFvalue, inputBarcode(end)];
                endValue = length(inputBarcode) - 1; 
            end
            for i = 2:endValue;
                if inputBarcode(i) == 0 && inputBarcode(i-1) > 0;
                   binsToSetToPF = [binsToSetToPF, i];
                   PFvalue = [PFvalue, inputBarcode(i-1)];
                end
            end
            inputBarcode(binsToSetToPF) = PFvalue; 
            % Repeat the same thing again to account for the next bin in. 
            endValue = length(inputBarcode);
            if inputBarcode(end) > 0 && inputBarcode(1) == 0; 
                binsToSetToPF = [binsToSetToPF, 1]; 
                PFvalue = [PFvalue, inputBarcode(end)];
                endValue = length(inputBarcode) - 1; 
            end
            for i = 2:endValue; 
                if inputBarcode(i) == 0 && inputBarcode(i-1) > 0;
                   binsToSetToPF = [binsToSetToPF, i];
                   PFvalue = [PFvalue, inputBarcode(i-1)];
                end
            end
            inputBarcode(binsToSetToPF) = PFvalue; 

            % If the field now 'straddles' the track
            % (meaning, if it is goes past the beginning or
            % the end) shift the track so that this isn't
            % the case
            if inputBarcode(1) ~=0 && inputBarcode(end) ~=0; 
                zeroIndices = find(inputBarcode == 0); 
                inputBarcode = circshift(inputBarcode, -zeroIndices(1), 2); 
                barcodeIndices = circshift(barcodeIndices, -zeroIndices(1), 2);
            end

            % Some of the fields may have merged, account for this by running the
            % bwlabel one more time
            newBarcode = bwlabel(inputBarcode); 
            if max(newBarcode) ~= max(inputBarcode); 
                inputBarcode = newBarcode; 
            end
            
        elseif modification < 0; % make the field smaller
            for iNumBins = 1:abs(modification);
                for i = 1:length(inputBarcode)-1; 
                   if inputBarcode(i) == 0 && inputBarcode(i+1) > 0; 
                       binsToSetToPF = [binsToSetToPF, i+1];
                   end
                end
                % If the current bin is 0 and the last one was > 0, set the last bin to 0
                for i = 2:length(inputBarcode); 
                    if inputBarcode(i) == 0 && inputBarcode(i-1) > 0;
                       binsToSetToPF = [binsToSetToPF, i-1];
                    end
                end
                inputBarcode(binsToSetToPF) = 0; 
            end
        end
        
        outputBarcode = inputBarcode;
        
    end

end
