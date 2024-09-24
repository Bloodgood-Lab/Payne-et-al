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
    for iGenotype = 1%:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype});
        for iFR = 1%:length(fieldnames(genotypeData)); 
            FRoptions = fieldnames(genotypeData); 
            FRdata = genotypeData.(FRoptions{iFR}); 
            for iAnimal = 1%:length(FRdata); 
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
                            for iDir = 1%:length(directions);
                                % Extract data based on running direction
                                outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir), 'spatialMetrics');
                                map = outputData.map; 
                                
                                % Assign data based on running direction
                                if strcmp(directions(iDir), 'cw') == 1; 
                                    [barcode, PFsize, PFnumber] = getPlaceFields(map, thresholds);
                                    % For use in later control analyses, modify the 
                                    % barcode to make the fields larger or smaller
                                    [biggerFieldBarcode, biggerFieldIndices] = modifyPlaceFields(barcode, biggerFieldModification);
                                    [smallerFieldBarcode, smallerFieldIndices] = modifyPlaceFields(barcode, smallerFieldModification);
                                    % Define data within the structure
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFsize.cw = PFsize;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFnumber.cw = PFnumber;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.original.cw = barcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.barcode.cw = biggerFieldBarcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.indices.cw = biggerFieldIndices;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.barcode.cw = smallerFieldBarcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.indices.cw = smallerFieldIndices;
                                elseif strcmp(directions(iDir), 'ccw') == 1; 
                                    [barcode, PFsize, PFnumber] = getPlaceFields(map, thresholds);
                                    % For use in later control analyses, modify the 
                                    % barcode to make the fields larger or smaller
                                    [biggerFieldBarcode, biggerFieldIndices] = modifyPlaceFields(barcode, biggerFieldModification);
                                    [smallerFieldBarcode, smallerFieldIndices] = modifyPlaceFields(barcode, smallerFieldModification);
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFsize.ccw = PFsize;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFnumber.ccw = PFnumber;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.ccw = barcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.barcode.ccw = biggerFieldBarcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.indices.ccw = biggerFieldIndices;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.barcode.ccw = smallerFieldBarcode;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.bigger.indices.ccw = smallerFieldIndices;
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
    
    function [barcode, PFSize, PFNumber] = getPlaceFields(map, thresholds)
        % Identifies the bins that are in-field
        % Inputs: 
        %   1) map: the collapsed map to be analyzed
        %   2) settings: the settings file which should contain lowThresh: the
        %      lower threshold for what is considered a place field and
        %      highThresh: the maximum firing rate that the field should
        %      contain
        % Output: 
        %   1) barcode: the barcode that tells you which bins belong to each
        %      field
        %   2) PFSize: the size of place fields
        %   3) PFNumber: the number of place fields
        % Steps: 
        %   1) Find the bins above the low threshold 
        %   2) Only include place fields if the max is above the high threshold
        %   3) Re-number the fields so that they are ordered by firing rate 
        %      (i.e. field 1 has the highest FR, then field 2, etc.)
        %   4) Get the number of place fields 

        %% Step 1: Find the bins above the low threshold 
        low_threshold = thresholds(1) * nanmax(map(:)); 
        barcode = bwlabel(map > low_threshold);

        %% Step 2: Only include fields where the max is above the high thresh
        high_threshold = thresholds(2) * nanmax(map(:));
        tempSize = []; FR = []; count = 1; 
        for i = 1:max(max(barcode)); 
            if max(map(barcode == i)) >= high_threshold;
                tempSize(count) = 4*length(find(barcode == i)); % Multiply by 4 since bin size is 4 cm.
                FR(count) = 16000 * nanmax(map(barcode == i)); % Get the max firing rate so that the fields can be ordered by FR
                barcode(barcode == i) = count; 
                count = count + 1;
            elseif max(map(barcode == i)) < high_threshold;  
                barcode(barcode == i) = 0; 
            end;
        end

        %% Step 3: Reorder the fields in order of firing rate
        [~, sortedInd] = sort(FR, 'descend'); % Get the order of the fields from highest FR to lowest
        barcode = 10*barcode; % Need this to be a larger value so that you can reorganize and reset values starting at 1
        PFSize = []; 
        for iField = 1:length(sortedInd); 
            barcode(barcode == 10*sortedInd(iField)) = iField;
            PFSize(iField) = tempSize(sortedInd(iField));
        end

        %% Step 4: Get the number of place fields
        PFNumber = length(PFSize);
    end

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
