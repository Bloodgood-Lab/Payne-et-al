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
                                outputData = assignVariableByDirection_v1_20240905(FRdata{iAnimal}(iCluster), directions(iDir));
                                map = outputData.map; 
                                
                                % Assign data based on running direction
                                if strcmp(directions(iDir), 'cw') == 1; 
                                    [barcode, PFsize, PFnumber] = getPlaceFields(map, thresholds);
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFsize.cw = PFsize;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFnumber.cw = PFnumber;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.cw = barcode;
                                elseif strcmp(directions(iDir), 'ccw') == 1; 
                                    [barcode, PFsize, PFnumber] = getPlaceFields(map, thresholds);
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFsize.ccw = PFsize;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.PFnumber.ccw = PFnumber;
                                    data.(genotypes{iGenotype}).(FRoptions{iFR}){iAnimal}(iCluster).spatialMetrics.barcode.ccw = barcode;
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
end
