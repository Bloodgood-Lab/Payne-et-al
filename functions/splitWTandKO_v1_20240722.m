function filePaths = splitWTandKO_v1_20240722(mainDataStructure) 
    % Splits the data into WT and KO
    % Written by Anja Payne
    % Last Modified: 07/22/2024

    % Inputs:
    %   1) mainDataStructure: the matlab structure where the file paths are
    %      stored

    % Outputs:
    %   1) filePaths: A structure split into WT and KO that contains the 
    %      file paths for the cells

    % Steps:
    %   1) Split into WT and KO based off of their opto response ('yes' or
    %      'no')
    %   2) Give the user the chance to save the data

    
    %% Step 1: Split into WT and KO based off of their opto response 
    filePaths = {}; WTcount = 1; KOcount = 1; 
    for iAnimal = 1:length(mainDataStructure); 
        for iCell = 1:length(mainDataStructure{iAnimal});
            if strcmp(mainDataStructure{iAnimal}(iCell).opto, 'No') == 1; 
                filePaths.WT(WTcount) = mainDataStructure{iAnimal}(iCell); 
                WTcount = WTcount + 1;
            end
            if strcmp(mainDataStructure{iAnimal}(iCell).opto, 'Yes') == 1; 
                filePaths.KO(KOcount) = mainDataStructure{iAnimal}(iCell); 
                KOcount = KOcount + 1;
            end
        end
    end
    
    %% Step 2: Give the user the chance to save the data
    saveFile_v1_20240718(filePaths, 'filePathsByGenotype') 
    
