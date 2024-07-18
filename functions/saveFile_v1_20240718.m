function saveFile_v1_20240718(data, dataName) 
    % Given user-defined save folder, save the data 
    % Written by Anja Payne
    % Last Modified: 07/18/2024

    % Inputs:
    %   1) saveData: data to save
    %   2) dataName: name of data to save
    %   2) also takes user-defined input from pop-up box

    % Outputs: None

    % Steps:
    %   1) Ask the user if they want to save the data
    %   2) If they say yes, have them select where
    %   3) If the file already exists in a different version, create a new
    %      version with today's date
    
    %% Step 1: Ask the user if they want to save the data
    saveChoice = questdlg('Do you want to save the data?', ...
    'Choose Option', ...
    'Save Data', 'Cancel', 'Cancel');

    switch saveChoice
        %% Step 2: If they say yes, have them select where
        case 'Save Data'
            fileNameBase = dataName; 
            savePathName = uigetdir();
            currentFiles = dir(savePathName);
            currentFiles = currentFiles(~ismember({currentFiles.name}, {'.', '..'}));
            [m, ~] = size(currentFiles); 
            if m == 0; 
                saveFile = [fileNameBase, '_v01_', date, '.mat'];
                save([savePathName, '\', saveFile], 'data'); 
            %% Step 3: If the file already exists, create a new version
            elseif m > 0;
                for iFiles = 1:m; 
                    % Check to see if the file already exists
                    fileNameBase_toCompare = currentFiles(iFiles).name(1:length(fileNameBase));
                    if strcmp(fileNameBase_toCompare, fileNameBase) == 1
                        oldVersion = currentFiles(iFiles).name(length(fileNameBase)+3:length(fileNameBase)+4);
                        newVersion = str2double(oldVersion) + 1;
                        newVersion = sprintf('%02d', newVersion);
                        saveFile = [fileNameBase, '_v', num2str(newVersion), '_', date, '.mat'];
                        save([savePathName, '\', saveFile], 'data'); 
                    end
                end
            end
        case 'Cancel'
            disp('Operation canceled by the user'); 
        otherwise
            error('Unexpected option selected.');   
    end