function [barcode, PFSize, PFNumber] = getPlaceFields_v1_20250425(map, thresholds)
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