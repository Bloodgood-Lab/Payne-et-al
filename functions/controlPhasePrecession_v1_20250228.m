function data = controlPhasePrecession_v1_20250228(data, settings)
    % Generates plots related to phase precession analysis
    % Written by Anja Payne
    % Last Modified: 03/11/2025
    
    % Inputs:
    %   1) data: the matlab structure where the in-field theta phases
    %      by-trial are saved
    %   2) settings: the settings file where the settings for the phase
    %      precession calculations are saved

    % Outputs:
    %   
    %     
    
    % Steps/Figures:
    %   1) Get the average size of the field for each cell that was
    %      included in the phase precession slope analysis
    %   2) Match the distributions of the place field sizes then compare 
    %      the slopes
    
    
    %% Step 1: Get the average place field size
    for iGenotype = 1:length(fieldnames(data.cellData));
        genotypes = fieldnames(data.cellData); 
        genotypeData = data.cellData.(genotypes{iGenotype}); 
        allFieldSizes = []; 
        
        % Run analysis for high-firing cells only
        FRdata = genotypeData.highFiring;
        for iAnimal = 1:length(FRdata); 
            % Skip if empty
            if isempty(FRdata{iAnimal}) == 1; 
                continue
            else
                [~,n] = size(FRdata{iAnimal});
                for iCluster = 1:n;
                    % Skip if empty
                    if isempty(FRdata{iAnimal}(iCluster).metaData) == 1; 
                        display(['Cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal), ' is empty, skipping']);
                        continue
                    else
                        display(['Calculating for cluster ', num2str(iCluster) ' of animal ', num2str(iAnimal)]);
                        
                        % Assign variables based on running direction
                        directions = fieldnames(FRdata{iAnimal}(iCluster).spatialMetrics.barcode.original);
                        for iDir = 1:length(directions);
                            % Loop through fields
                            fieldsToAnalyze = settings.phasePrecession.fieldsToAnalyze; 
                            if strcmp(fieldsToAnalyze, 'all fields') == 1;
                                numField = length(FRdata{iAnimal}(iCluster).spatialMetrics.spatialMetrics.PFsize.(directions{iDir})); 
                            elseif strcmp(fieldsToAnalyze, 'best field') == 1;
                                numField = 1; 
                            end
                            
                            for iField = 1:numField;
                                singleFieldSize = FRdata{iAnimal}(iCluster).spatialMetrics.PFsize.(directions{iDir})(iField);  
                                allFieldSizes = [allFieldSizes, singleFieldSize];
                            end
                        end
                    end
                end
            end
        end
        data.populationData(iGenotype).phasePrecession.averageFieldSizes = allFieldSizes;
    end
    
    
    
    %% Step 2: Match the distributions of the place field sizes then compare the slopes
    
    
    
    
    
    
end
