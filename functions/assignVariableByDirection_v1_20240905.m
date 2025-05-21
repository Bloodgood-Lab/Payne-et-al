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
        if strcmp(analysisType, 'spikeTiming') == 1; 
            outputData.spkTimes = inputData.highVelocityData.spikeTimesByTrial.cw;
        end
        if strcmp(analysisType, 'spatialMetrics') == 1;
            outputData.map = inputData.rateMaps.trialAverageMap.cw;
            outputData.fullMap = inputData.rateMaps.rateMap.cw; 
            outputData.timeMap = inputData.rateMaps.timeMap.cw; 
        end
        if strcmp(analysisType, 'getInField') == 1;
            outputData.barcode = inputData.spatialMetrics.barcode.original.cw;
            outputData.binnedSpkPosForInField = inputData.binnedSpikesByTrial.highVelocityData.binnedSpkPos.cw;
            outputData.spkPosForInField = inputData.binnedSpikesByTrial.highVelocityData.unbinnedSpkPos.cw;
            outputData.spkTimes = inputData.binnedSpikesByTrial.highVelocityData.binnedSpkTimes.cw;
            outputData.bursts = inputData.spikeTiming.bursts.cw; 
            outputData.singles = inputData.spikeTiming.singles.cw; 
        end
        if strcmp(analysisType, 'theta') == 1;
            outputData.spikesByDirection = inputData.inField.inFieldSpkTimes.cw; 
            outputData.inFieldBursts = inputData.inField.inFieldBursts.cw; 
            outputData.inFieldSingles = inputData.inField.inFieldSingles.cw; 
        end
        if strcmp(analysisType, 'phasePrecession') == 1;
            outputData.spikesByDirection = inputData.inField.inFieldSpkTimes.cw; 
            outputData.spkPhs = inputData.theta.allSpikes.phases.cw; 
            outputData.spkPos = inputData.inField.inFieldSpkPos.cw; 
            outputData.binnedSpkPos = inputData.inField.inFieldBinnedSpkPos.cw; 
        end
        if strcmp(analysisType, 'plotPhasePrecession') == 1; 
            outputData.spikesByDirection = inputData.inField.inFieldSpkTimes.cw;
            outputData.binnedSpkPos = inputData.inField.inFieldBinnedSpkPos.cw; 
            outputData.spkPhsForPlot = inputData.phasePrecession.phsInput.cw;
            outputData.spkPosForPlot = inputData.phasePrecession.posInput.cw;
            outputData.slopeMedian = inputData.phasePrecession.medianSlope.cw;
            outputData.allSlopes = inputData.phasePrecession.allSlopes.cw;
            outputData.offsets = inputData.phasePrecession.offsets.cw; 
            outputData.rSquared = inputData.phasePrecession.rSquared.cw; 
            outputData.pValues = inputData.phasePrecession.pValues.cw; 
        end
        
    elseif strcmp(dir, 'ccw') == 1; 
        if strcmp(analysisType, 'spikeTiming') == 1; 
            outputData.spkTimes = inputData.highVelocityData.spikeTimesByTrial.ccw;
        end
        if strcmp(analysisType, 'spatialMetrics') == 1;
            outputData.map = inputData.rateMaps.trialAverageMap.ccw;
            outputData.fullMap = inputData.rateMaps.rateMap.ccw; 
            outputData.timeMap = inputData.rateMaps.timeMap.ccw; 
        end
        if strcmp(analysisType, 'getInField') == 1;
            outputData.barcode = inputData.spatialMetrics.barcode.original.ccw;
            outputData.binnedSpkPosForInField = inputData.binnedSpikesByTrial.highVelocityData.binnedSpkPos.ccw;
            outputData.spkPosForInField = inputData.binnedSpikesByTrial.highVelocityData.unbinnedSpkPos.ccw;
            outputData.spkTimes = inputData.binnedSpikesByTrial.highVelocityData.binnedSpkTimes.ccw;
            outputData.bursts = inputData.spikeTiming.bursts.ccw; 
            outputData.singles = inputData.spikeTiming.singles.ccw;         
        end
        if strcmp(analysisType, 'theta') == 1;
            outputData.spikesByDirection = inputData.inField.inFieldSpkTimes.ccw; 
            outputData.inFieldBursts = inputData.inField.inFieldBursts.ccw; 
            outputData.inFieldSingles = inputData.inField.inFieldSingles.ccw;
        end
        if strcmp(analysisType, 'phasePrecession') == 1;
            outputData.spikesByDirection = inputData.inField.inFieldSpkTimes.ccw; 
            outputData.spkPhs = inputData.theta.allSpikes.phases.ccw; 
            outputData.spkPos = inputData.inField.inFieldSpkPos.ccw; 
            outputData.binnedSpkPos = inputData.inField.inFieldBinnedSpkPos.ccw;
        end
        if strcmp(analysisType, 'plotPhasePrecession') == 1; 
            outputData.spikesByDirection = inputData.inField.inFieldSpkTimes.ccw;
            outputData.binnedSpkPos = inputData.inField.inFieldBinnedSpkPos.ccw;
            outputData.spkPhsForPlot = inputData.phasePrecession.phsInput.ccw;
            outputData.spkPosForPlot = inputData.phasePrecession.posInput.ccw;
            outputData.slopeMedian = inputData.phasePrecession.medianSlope.ccw;
            outputData.allSlopes = inputData.phasePrecession.allSlopes.ccw;
            outputData.offsets = inputData.phasePrecession.offsets.ccw; 
            outputData.rSquared = inputData.phasePrecession.rSquared.ccw; 
            outputData.pValues = inputData.phasePrecession.pValues.ccw; 
        end
    end
    
    
                                    