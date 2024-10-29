function outputData = calculatePhasePrecession_v1_20241020(inputData, settings)
    % Calculates the phase precession depending on the mode
    % Written by Anja Payne
    % Last Modified: 10/20/2024
    
    % Inputs:
    %   1) inputData: the data for that specific field
    %   2) settings: the settings file
    
    % Outputs:
    %   1) outputData: the extracted variables for that direction
    
   
    if strcmp(settings.phasePrecession.circularity, 'shift') == 1;
        % To account for circularity in the data, find the phase
        % shift that accounts for the best correlation between
        % position and spike theta phases
        spkPhsDegrees = rad2deg(inputData.spkPhs);
        phasesToShiftBy = [0:180/36:360];
        correlation = [];
        for iShift = 1:length(phasesToShiftBy);
            testShiftedPhs = circshift(spkPhsDegrees, phasesToShiftBy(iShift));
            R = corrcoef(inputData.spkPos, testShiftedPhs); correlation(iShift) = R(2);
        end
        [~, maxInd] = min(correlation);
        shiftedPhs = circshift(spkPhsDegrees, phasesToShiftBy(maxInd));
        outputData.spkPhs = deg2rad(shiftedPhs) + pi;
        outputData.spkPos = inputData.spkPos; 
    elseif strcmp(settings.phasePrecession.circularity, 'double') == 1;
        outputData.spkPhs = [inputData.spkPhs+pi; inputData.spkPhs+3*pi];
        outputData.spkPos = [inputData.spkPos; inputData.spkPos];
    elseif strcmp(settings.phasePrecession.circularity, 'none') == 1;
        outputData.spkPos = inputData.spkPos;
        outputData.spkPhs = inputData.spkPhs+pi;
        R = corrcoef(outputData.spkPos, outputData.spkPhs); correlation = R(2); maxInd = 1; 
    end

    % Get the phase precession for that trial
    [cir, lin, circFit] = thetaPrecess(outputData.spkPhs, outputData.spkPos, settings.phasePrecession.slopeRange); 
    if strcmp(settings.phasePrecession.fit, 'circularSlope') == 1;
        outputData.trialSlope = cir.Alpha; 
        outputData.trialOffset = cir.Phi0;
        outputData.r2 = (cir.Coeff)^2;
        outputData.p = cir.pValue;
    elseif strcmp(settings.phasePrecession.fit, 'linearFit') == 1;
        outputData.trialSlope = lin.Alpha;
        outputData.trialOffset = lin.Phi0;
        outputData.r2 = (lin.r)^2;
        outputData.p = lin.p;
    elseif strcmp(settings.phasePrecession.fit, 'circularFit') == 1;
        outputData.trialSlope = circFit.Alpha;
        outputData.trialOffset = circFit.Phi0;
        outputData.r2 = (circFit.rho)^2;   
        outputData.p = circFit.p; 
    end
                                    