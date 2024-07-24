function data = getSpatialMetrics_v1_20240724(data, settings)
    % Gets the spatial metrics from the trial-averaged rate maps
    % Written by Anja Payne
    % Last Modified: 07/24/2024

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
    %      average size, and number of place fields
    %   2) Ask the user if they want to save the newly generated data
    
    
    %% Step 1: Get the barcode, average size, and number of place fields
    thresholds(1) = settings.rateMap.lowThresh; 
    thresholds(2) = settings.rateMap.highThresh;
    for iGenotype = 1:length(fieldnames(data));
        genotypes = fieldnames(data); 
        genotypeData = data.(genotypes{iGenotype}); 
        for iAnimal = 1:length(genotypeData); 
            if isempty(genotypeData{iAnimal}) == 1; 
                continue
            else
                for iCluster = 1:length(genotypeData{iAnimal}.rateMaps);
                    % Direction 1
                    map = genotypeData{iAnimal}.rateMaps{iCluster}.dir1;
                    [barcode, PFsize, PFnumber] = getPlaceFields(map, thresholds);
                    data.(genotypes{iGenotype}){iAnimal}.PFsize{iCluster}.dir1 = PFsize;
                    data.(genotypes{iGenotype}){iAnimal}.PFnumber{iCluster}.dir1 = PFnumber;
                    data.(genotypes{iGenotype}){iAnimal}.barcode{iCluster}.dir1 = barcode;
                    
                    % Direction 2
                    map = genotypeData{iAnimal}.rateMaps{iCluster}.dir2;
                    [barcode, PFsize, PFnumber] = getPlaceFields(map, thresholds);
                    data.(genotypes{iGenotype}){iAnimal}.PFsize{iCluster}.dir2 = PFsize;
                    data.(genotypes{iGenotype}){iAnimal}.PFnumber{iCluster}.dir2 = PFnumber;
                    data.(genotypes{iGenotype}){iAnimal}.barcode{iCluster}.dir2 = barcode;
                end
            end
        end
    end
    
    %% Step 2: 
    saveFile_v1_20240718(data, settings, 'spatialMetrics') 

    
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
end
