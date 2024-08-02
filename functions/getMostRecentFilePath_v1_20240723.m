function filePath = getMostRecentFilePath_v1_20240723(fileNameBase) 
    % Given user-defined folder and base file name, get the most recent
    % version
    % Written by Anja Payne
    % Last Modified: 07/31/2024

    % Inputs:
    %   1) fileNameBase: file base name
    %   2) also takes user-defined folder location from pop-up box

    % Outputs: 
    %   1) filePath{1}: the user-defined pathway to where the processed
    %      data is stored
    %   2) filePath{2}: the most recent version, saved as a string
    %   3) filePath{3}: any appended information (usually the date)
    %      associated with the most recent version

    % Steps:
    %   1) Get the folder from the user
    %   2) Find the most recent file with the specified base name in that
    %      folder
    
    %% Step 1: Get the folder from the user
    userSelectedFolder = uigetdir();
    filePath{1} = [userSelectedFolder, '\', fileNameBase];
    
    %% Step 2: Find the most recent file with the specified base name 
    currentFiles = dir(filePath{1});
    currentFiles = currentFiles(~ismember({currentFiles.name}, {'.', '..'}));
    [m, ~] = size(currentFiles);
    if m == 0; 
    	filePath{2} = '00'; 
        filePath{3} = date; 
    elseif m > 0;
        allVersions = []; allAppended = {};
        for iFiles = 1:m; 
            % If the file is shorter than the fileNameBase, skip it
            if length(currentFiles(iFiles).name) < length(fileNameBase);
                continue;
            else
                % Check to see if the file already exists
                fileNameBase_toCompare = currentFiles(iFiles).name(1:length(fileNameBase));
                if strcmp(fileNameBase_toCompare, fileNameBase) == 1
                    version = currentFiles(iFiles).name(length(fileNameBase)+3:length(fileNameBase)+4);
                    append = currentFiles(iFiles).name(length(fileNameBase)+5:end);
                    allVersions(iFiles) = str2double(version);
                    allAppended{iFiles} = append; 
                end
            end
        end
        if isempty(allAppended) == 1; 
            filePath{2} = '00'; 
            filePath{3} = ['_', date]; 
        else
            [maxVersion, maxI] = nanmax(allVersions);
            maxVersion = sprintf('%02d', maxVersion); 
            maxAppended = allAppended{maxI}; 
            filePath{2} = maxVersion;
            filePath{3} = maxAppended;
        end
    end
 
    
        
