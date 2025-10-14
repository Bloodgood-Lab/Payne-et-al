function [map, timeMap, spikeMap] = calculateLinearizedRateMap(pos, spikes, nBins)
    % Gets the rate map for a linearized representation of the map
    % Original code from Leutgeb lab with modifications by A. Payne
    
    % Inputs:
    %   1) pos: the binned positions
    %   2) spikes: the binned spikes
    % Outputs:
    %   1) map: a heatmap of firing rate per bin
    %   2) timeMap: a heatmap of time per bin
    %   3) spikeMap: a heatmap of spike number per bin
    
    spikeMap = zeros(1, nBins); 
    for iSpike = 1:length(spikes);
        spikeMap(spikes(iSpike)) = spikeMap(spikes(iSpike)) + 1; 
    end
    timeMap = zeros(1, nBins); 
    for iTime = 1:length(pos);
       timeMap(pos(iTime)) = timeMap(pos(iTime)) + 1; 
    end

    % Set unoccupied pixels to NaN; 
    minOccupancy = 1; 
    timeMap(timeMap < minOccupancy) = NaN; 
    rawMap = spikeMap./timeMap; 
    % Shift the map so that the reward zone is at the end
    rawMap = circshift(rawMap, -2, 2);
    % Add some padding for smoothing to work effectively
    rawMap = padarray(rawMap', 2, NaN); 
    
    % Smooth the rate map using the Gaussian filter as well as an adaptive
    % smoothing technique
    % Define the Gaussian smoothing that will be applied to each bin. 
    box = [0.0200 0.1000 0.1600 0.1000 0.0200]; 
    for i = 3:length(rawMap)-2; 
        current_sum = 0;
        current_box = 0;
        box_i = 0;
        for ii = i-2:i+2;
          box_i = box_i + 1;
          if ~isnan(rawMap(ii));
             current_sum = current_sum + rawMap(ii) * box(box_i);
             current_box = current_box + box(box_i); 
          end
        end
        map(i) = current_sum/current_box; 
    end
    
    map = map(3:end); 
    
end
