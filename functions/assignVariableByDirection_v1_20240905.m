function outputData = assignVariableByDirection_v1_20240905(inputData, dir, analysisType)
    % Given data from one cell obtained from the data structure, extract
    % the appropriate data based on the direction
    % Written by Anja Payne
    % Last Modified: 09/05/2024

    % Inputs:
    %   1) inputData: the data for that specific cell
    %   2) dir: either direction 1 or 2
    
    % Outputs:
    %   1) outputData: the extracted variables for that direction

    if strcmp(dir, 'cw') == 1; 
        outputData.map = inputData.rateMaps.trialAverageMap.cw;
        outputData.barcode = inputData.spatialMetrics.barcode.cw;
        outputData.binnedSpkPosForInField = inputData.binnedSpikesByTrial.allVelocities.binnedSpkPos.cw;
        outputData.spkPosForInField = inputData.binnedSpikesByTrial.allVelocities.unbinnedSpkPos.cw;
        outputData.spkTimes = inputData.binnedSpikesByTrial.allVelocities.binnedSpkTimes.cw;
        outputData.spikesByDirection = inputData.inField.inFieldSpkTimes.cw; 
        outputData.spkPhs = inputData.theta.phases.cw; 
        outputData.spkPos = inputData.inField.inFieldSpkPos.cw; 
        outputData.binnedSpkPos = inputData.inField.inFieldBinnedSpkPos.cw; 
        if strcmp(analysisType, 'plotPhasePrecession') == 1; 
            outputData.spkPhsForPlot = inputData.phasePrecession.phsInput.cw;
            outputData.spkPosForPlot = inputData.phasePrecession.posInput.cw;
            outputData.lineFit = inputData.phasePrecession.fitInfo.cw;
            outputData.slopeMedian = inputData.phasePrecession.medianSlope.cw;
            outputData.allSlopes = inputData.phasePrecession.allSlopes.cw;
        end
        
    elseif strcmp(dir, 'ccw') == 1; 
        outputData.map = inputData.rateMaps.trialAverageMap.ccw;
        outputData.barcode = inputData.spatialMetrics.barcode.ccw;
        outputData.binnedSpkPosForInField = inputData.binnedSpikesByTrial.allVelocities.binnedSpkPos.ccw;
        outputData.spkPosForInField = inputData.binnedSpikesByTrial.allVelocities.unbinnedSpkPos.ccw;
        outputData.spkTimes = inputData.binnedSpikesByTrial.allVelocities.binnedSpkTimes.ccw;
        outputData.spikesByDirection = inputData.inField.inFieldSpkTimes.ccw; 
        outputData.spkPhs = inputData.theta.phases.ccw; 
        outputData.spkPos = inputData.inField.inFieldSpkPos.ccw; 
        outputData.binnedSpkPos = inputData.inField.inFieldBinnedSpkPos.ccw;
        if strcmp(analysisType, 'plotPhasePrecession') == 1; 
            outputData.spkPhsForPlot = inputData.phasePrecession.phsInput.ccw;
            outputData.spkPosForPlot = inputData.phasePrecession.posInput.ccw;
            outputData.lineFit = inputData.phasePrecession.fitInfo.ccw;
            outputData.slopeMedian = inputData.phasePrecession.medianSlope.ccw;
            outputData.allSlopes = inputData.phasePrecession.allSlopes.ccw;
        end
    end
    
    
                                    