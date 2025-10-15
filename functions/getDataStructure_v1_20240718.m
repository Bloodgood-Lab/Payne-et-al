function mainDataStructure = getDataStructure_v1_20240718(excelFiles, excelFolder) 
    % Get the data structure by pulling relevant information from excel 
    % Written by Anja Payne
    % Last Modified: 07/18/2024

    % Inputs:
    %   1) excelFiles: the list of excel files where the data is stored
    %   2) excelFolder: the path to the folder where the excel files are stored

    % Outputs:

    % Steps:
    %   1) Ask the user if they want to load the data or recalculate it
    %   2) If the user wants to recalculate the data, convert the excel 
    %      file to matlab file. Otherwise load data.
    %   3) Ask the user if they want to save the newly generated data

    
    %% Step 1: Ask the user if they want to load the data or recalculate it

    loadChoice = questdlg('Do you want to load a pre-existing data structure or calculate a new one?', ...
        'Choose Option', ...
        'Load Data', 'Calculate New Structure', 'Cancel', 'Cancel');
   
    %% Step 2: Either load or recalculate based on user input
    switch loadChoice
        case 'Load Data' 
            [filename, pathname] = uigetfile();
            load([pathname, filename]);
            mainDataStructure = data; 
        case 'Calculate New Structure'
            mainDataStructure = {}; 
            for iCohort = 1:length(excelFiles); 
                excelFilePath = [excelFolder, excelFiles{iCohort}, '.xlsx']; 

                if exist(excelFilePath, 'file') ~=2;
                    error('File "%s" does not exist.', excelFilePath); 
                end

                cohortDataStructure = convertExcelToMat(excelFilePath);
                mainDataStructure = [mainDataStructure, cohortDataStructure]; 
            end
            % Create an empty settings structure
            settings = struct();
            
            %% Step 3: Ask the user if they want to save the newly generated data
            saveFile_v1_20240718(mainDataStructure, settings, 'mainDataStructure') 
        
        case 'Cancel'
            disp('Operation canceled by the user'); 
        otherwise
            error('Unexpected option selected.');
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function animalInfo = convertExcelToMat(excelFile)
    % Loads an excel workbook into matlab
    % Inputs:
    %   1) excelFile: the path to the excel file
    % Outputs:
    %   1) animalInfo: the data structure
    % Steps:
    %   1) Load the sheet names so it knows which animals to iterate over (each
    %      sheet is one animal)
    %   2) If the session is 'track' type, save the 'useful' data (in our case 
    %      this is the directory and the opto-tagging result) as a structure 
    %      organized according to structure{iAnimal}{iCell}{filepath, optotag result}

    %% Step 1:
    [~, animals] = xlsfinfo(excelFile);

    %% Step 2:
    for iAnimal = 1:length(animals)
        [~, data] = xlsread(excelFile, iAnimal);
        [m, ~] = size(data); 
        cellCount = 1; 
        for iCluster = 2:m;
            if strcmp(data{iCluster, 2}, 'Track') == 1; 
                animalInfo{iAnimal}(cellCount).directory = data(iCluster, 1);
                animalInfo{iAnimal}(cellCount).fileName = data(iCluster, 3);
                animalInfo{iAnimal}(cellCount).opto = data(iCluster, 6);
                cellCount = cellCount + 1; % iterate
            end
        end
    end
end
