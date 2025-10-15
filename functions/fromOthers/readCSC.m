function [cscData, header] = readCSC(fileName)
% Function obtained from the Leutgeb Lab

% Set the field selection for reading CSC files. 1 = Add parameter, 0 = skip
% parameter
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 1; % Channel number
fieldSelection(3) = 1; % Sample Frequency
fieldSelection(4) = 1; % Number of valid samples
fieldSelection(5) = 1; % Samples
extractHeader = 1; %1 = Yes, 0 = No.
% 5 different extraction modes, see help file for Nlx2MatCSC_v3
extractMode = 1; % Extract all data

[ts,chNum,sampFreq,numSamp,samp, header] = Nlx2MatCSC(fileName,fieldSelection,extractHeader,extractMode);

samp = reshape(samp,[],1);
%disp(['Size of sample: ',num2str(length(samp))]) %display
% Put the data in a cell array
cscData = struct('ts',ts,'chNum',chNum,'sampFreq',sampFreq,...
       'numSamp',numSamp,'samp',samp);

% if combined_samples > 0
%     n_intervals = fix(size(trace.samp,2)/combined_samples);
%     EEG = reshape(trace.samp(1:512,1:combined_samples*n_intervals),combined_samples*512,n_intervals);
% else
%     
% end



%number of sampling points,number of sets
%end function
% Reshape into column vector

