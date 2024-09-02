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
        % Create the main folder if needed
        mainFolder = [settings.filePath{1}, '\', settings.appendedFolder.name];
        if ~exist(mainFolder, 'dir')
            % Create the folder if it does not exist
            mkdir(mainFolder);
        end

        % Create the subfolder with the
        % version and today's data
        saveVersion = str2double(settings.filePath{2}) + 1;                                        
        saveVersion = sprintf('%02d', saveVersion);
        subFolder = [mainFolder, '\v', saveVersion, '_', date];
        if ~exist(subFolder, 'dir')
            % Create the folder if it does not exist
            mkdir(subFolder);
        end

        % Get the file name and save
        saveFile = [settings.name.geno, '_Animal', num2str(settings.name.animal),...
            '_Cluster', num2str(settings.name.cluster), '_', settings.name.direction];

        for iFileType = 1:length(settings.fileTypes)
            saveas(fig, [subFolder, '\', saveFile], settings.fileTypes{iFileType});
        end
    end
    