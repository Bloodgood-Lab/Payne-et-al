function saveFile_v1_20240718(processedDataPath, data, settings, dataName) 
    % Given user-defined save folder, save the data 
    % Written by Anja Payne
    % Last Modified: 07/31/2024

    % Inputs:
    %   1) saveData: data to save
    %   2) dataName: name of data to save
    %   2) also takes user-defined input from pop-up box

    % Outputs: None

    % Steps:
    %   1) Ask the user if they want to save the data
    %   2) If they say yes, have them select where and, if the file already
    %      exists in an earlier version, create a new version with today's
    %      date
    
    %% Step 1: Ask the user if they want to save the data
    saveChoice = questdlg('Do you want to save the data and (if applicable) the associated figures?', ...
    'Choose Option', ...
    'Yes', 'No', 'Cancel', 'Cancel');

    switch saveChoice
        %% Step 2: If they say yes, have them select where
        case 'Yes'
            % Ask the user to write a short message about the data
            prompt = {'Please enter a commit message:'};
            dlgTitle = 'Data Commit Message';
            numLines = [1, 50];
            defaultAnswer = {'No commit message entered for this data'};
            message = inputdlg(prompt, dlgTitle, numLines, defaultAnswer);
            
            % Using the same pathway for loading data, get the most recent
            % file and save the new data
            fileNameBase = dataName;
            filePath = getMostRecentFilePath_v1_20240723(fileNameBase, processedDataPath);
            savePathName = filePath{1};
            saveVersion = str2double(filePath{2}) + 1;
            saveVersion = sprintf('%02d', saveVersion); 
            saveFile = [fileNameBase, '_v', saveVersion, '_', date];
            if ~exist(savePathName, 'dir')
            	% Create the folder if it does not exist
            	mkdir(savePathName);
            end
            save([savePathName, '\', saveFile, '.mat'], 'data', 'settings', 'message'); 
            
        case 'No'
            disp('Data not saved'); 
        case 'Cancel'
            disp('Operation canceled by the user'); 
        otherwise
            error('Unexpected option selected.');   
    end