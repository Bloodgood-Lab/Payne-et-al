function saveFigure_v1_20240902(fig, settings)
    % Saves plots at the location and in the formats defined by settings
    % Last Modified: 09/02/2024

    % Inputs:
    %   1) fig: the figure handle for the figure your want to save
    %   2) settings: the settings file where the settings for saving are
    %      stored

    % Outputs:
    %   1) plots saved in the folder specified by the user
    
    if strcmp(settings.appendedFolder.binary, 'yes') == 1;
        
        % Create the subfolder with the version and today's date
        saveVersion = str2double(settings.filePath{2}) + 1;                                        
        saveVersion = sprintf('%02d', saveVersion);
        subFolder = [settings.filePath{1}, '\', settings.appendedFolder.name, '_v', saveVersion, '_', date];
        if ~exist(subFolder, 'dir')
            % Create the folder if it does not exist
            mkdir(subFolder);
        end

        % Get the file name and save
        saveFile = settings.name; 
        for iFileType = 1:length(settings.fileTypes)
            saveas(fig, [subFolder, '\', saveFile], settings.fileTypes{iFileType});
        end
        
        %{
        if ~exist([subFolder, '\commitMessage.txt'], 'file')
            % Ask the user to write a short message about the data
            prompt = {'Please enter a commit message:'};
            dlgTitle = 'Data Commit Message';
            numLines = [1, 50];
            defaultAnswer = {'No commit message entered for this figure'};
            commitMessage = inputdlg(prompt, dlgTitle, numLines, defaultAnswer);

            % Save the commit message
            commitFile = fopen([subFolder, '\commitMessage.txt'], 'w');
            fprintf(commitFile, '%s', commitMessage{1}); 
            fclose(commitFile); 
        end
        %}
        
    elseif strcmp(settings.appendedFolder.binary, 'no') == 1;
        
        if ~exist(settings.filePath{1}, 'dir')
            % Create the folder if it does not exist
            mkdir(settings.filePath{1});
        end
        
        % Create the subfolder with the version and today's data
        saveVersion = str2double(settings.filePath{2}) + 1;                                        
        saveVersion = sprintf('%02d', saveVersion);
        
        % Get the file name and save
        saveFile = settings.name; 
        for iFileType = 1:length(settings.fileTypes)
            saveas(fig, [settings.filePath{1}, '\', saveFile, '_v', saveVersion, '_', date], settings.fileTypes{iFileType});
        end
    end
    