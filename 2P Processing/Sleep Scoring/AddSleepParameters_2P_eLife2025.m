function [] = AddSleepParameters_2P_nNOS(mergedDataFileIDs,RestingBaselines,baselineType)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
for a = 1:size(mergedDataFileIDs,1)
    mergedDataFileID = mergedDataFileIDs(a,:);
    disp(['Adding sleep scoring parameters to ' mergedDataFileID '... (' num2str(a) '/' num2str(size(mergedDataFileIDs,1)) ')']); disp(' ')
    [~,~,fileDate,~,~,vesselID] = GetFileInfo2_2P_nNOS(mergedDataFileID);
    strDay = ConvertDate_2P_nNOS(fileDate);
    load(mergedDataFileID)
    specDataFileID = [mergedDataFileID(1:end - 14) 'SpecData.mat'];
    load(specDataFileID)
    %% BLOCK PURPOSE: Create folder for the Neural data of each electrode
    % hippocampal delta
    hippDelta = MergedData.data.hippocampalNeural.deltaBandPower;
    hippBaselineDelta = RestingBaselines.(baselineType).hippocampalNeural.deltaBandPower.(strDay);
    hippDeltaNeuro = (hippDelta - hippBaselineDelta)/hippBaselineDelta; 
    % cortical delta
    cortDelta = MergedData.data.corticalNeural.deltaBandPower;
    cortBaselineDelta = RestingBaselines.(baselineType).corticalNeural.deltaBandPower.(strDay);
    cortDeltaNeuro = (cortDelta - cortBaselineDelta)/cortBaselineDelta;
    % hippocampal theta
    hippTheta = MergedData.data.hippocampalNeural.thetaBandPower;
    hippBaselineTheta = RestingBaselines.(baselineType).hippocampalNeural.thetaBandPower.(strDay);
    hippThetaNeuro = (hippTheta - hippBaselineTheta)/hippBaselineTheta; 
    % cortical theta
    cortTheta = MergedData.data.corticalNeural.thetaBandPower;
    cortBaselineTheta = RestingBaselines.(baselineType).corticalNeural.thetaBandPower.(strDay);
    cortThetaNeuro = (cortTheta - cortBaselineTheta)/cortBaselineTheta;
    % hippocampal alpha
    hippAlpha = MergedData.data.hippocampalNeural.alphaBandPower;
    hippBaselineAlpha = RestingBaselines.(baselineType).hippocampalNeural.alphaBandPower.(strDay);
    hippAlphaNeuro = (hippAlpha - hippBaselineAlpha)/hippBaselineAlpha;  
    % cortical alpha
    cortAlpha = MergedData.data.corticalNeural.alphaBandPower;
    cortBaselineAlpha = RestingBaselines.(baselineType).corticalNeural.alphaBandPower.(strDay);
    cortAlphaNeuro = (cortAlpha - cortBaselineAlpha)/cortBaselineAlpha;
    % hippocampal beta
    hippBeta = MergedData.data.hippocampalNeural.betaBandPower;
    hippBaselineBeta = RestingBaselines.(baselineType).hippocampalNeural.betaBandPower.(strDay);
    hippBetaNeuro = (hippBeta - hippBaselineBeta)/hippBaselineBeta;  
    % cortical beta
    cortBeta = MergedData.data.corticalNeural.betaBandPower;
    cortBaselineBeta = RestingBaselines.(baselineType).corticalNeural.betaBandPower.(strDay);
    cortBetaNeuro = (cortBeta - cortBaselineBeta)/cortBaselineBeta;
    % hippocampal gamma
    hippGamma = MergedData.data.hippocampalNeural.gammaBandPower;
    hippBaselineGamma = RestingBaselines.(baselineType).hippocampalNeural.gammaBandPower.(strDay);
    hippGammaNeuro = (hippGamma - hippBaselineGamma)/hippBaselineGamma;
    % cortical gamma
    cortGamma = MergedData.data.corticalNeural.gammaBandPower;
    cortBaselineGamma = RestingBaselines.(baselineType).corticalNeural.gammaBandPower.(strDay);
    cortGammaNeuro = (cortGamma - cortBaselineGamma)/cortBaselineGamma;
    % hippocampal MUA
    hippMUA = MergedData.data.hippocampalNeural.muaPower;
    hippBaselineMUA = RestingBaselines.(baselineType).hippocampalNeural.muaPower.(strDay);
    hippMUANeuro = (hippMUA - hippBaselineMUA)/hippBaselineMUA;
    % cortical MUA
    cortMUA = MergedData.data.corticalNeural.muaPower;
    cortBaselineMUA = RestingBaselines.(baselineType).corticalNeural.muaPower.(strDay);
    cortMUANeuro = (cortMUA - cortBaselineMUA)/cortBaselineMUA;
    % Divide the neural signals into five second bins and put them in a cell array
    hippTempDeltaStruct = cell(180,1);
    hippTempThetaStruct = cell(180,1);
    hippTempAlphaStruct = cell(180,1);
    hippTempBetaStruct = cell(180,1);
    hippTempGammaStruct = cell(180,1);
    hippTempMUAStruct = cell(180,1);
    cortTempDeltaStruct = cell(180,1);
    cortTempThetaStruct = cell(180,1);
    cortTempAlphaStruct = cell(180,1);
    cortTempBetaStruct = cell(180,1);
    cortTempGammaStruct = cell(180,1);
    cortTempMUAStruct = cell(180,1);
    % loop through all samples across the 15 minutes in 5 second bins (180 total)
    for b = 1:180
        if b == 1
            % hippocampal
            hippTempDeltaStruct(b,1) = {hippDeltaNeuro(b:150)};
            hippTempThetaStruct(b,1) = {hippThetaNeuro(b:150)};
            hippTempAlphaStruct(b,1) = {hippAlphaNeuro(b:150)};
            hippTempBetaStruct(b,1) = {hippBetaNeuro(b:150)};
            hippTempGammaStruct(b,1) = {hippGammaNeuro(b:150)};
            hippTempMUAStruct(b,1) = {hippMUANeuro(b:150)};
            % cortical
            cortTempDeltaStruct(b,1) = {cortDeltaNeuro(b:150)};
            cortTempThetaStruct(b,1) = {cortThetaNeuro(b:150)};
            cortTempAlphaStruct(b,1) = {cortAlphaNeuro(b:150)};
            cortTempBetaStruct(b,1) = {cortBetaNeuro(b:150)};
            cortTempGammaStruct(b,1) = {cortGammaNeuro(b:150)};
            cortTempMUAStruct(b,1) = {cortMUANeuro(b:150)};
        elseif b == 180
            % hippocampal
            hippTempDeltaStruct(b,1) = {hippDeltaNeuro((((150*(b - 1)) + 1)):end)};
            hippTempThetaStruct(b,1) = {hippThetaNeuro((((150*(b - 1)) + 1)):end)};
            hippTempAlphaStruct(b,1) = {hippAlphaNeuro((((150*(b - 1)) + 1)):end)};
            hippTempBetaStruct(b,1) = {hippBetaNeuro((((150*(b - 1)) + 1)):end)};
            hippTempGammaStruct(b,1) = {hippGammaNeuro((((150*(b - 1)) + 1)):end)};
            hippTempMUAStruct(b,1) = {hippMUANeuro((((150*(b - 1)) + 1)):end)};
            % cortical
            cortTempDeltaStruct(b,1) = {cortDeltaNeuro((((150*(b - 1)) + 1)):end)};
            cortTempThetaStruct(b,1) = {cortThetaNeuro((((150*(b - 1)) + 1)):end)};
            cortTempAlphaStruct(b,1) = {cortAlphaNeuro((((150*(b - 1)) + 1)):end)};
            cortTempBetaStruct(b,1) = {cortBetaNeuro((((150*(b - 1)) + 1)):end)};
            cortTempGammaStruct(b,1) = {cortGammaNeuro((((150*(b - 1)) + 1)):end)};
            cortTempMUAStruct(b,1) = {cortMUANeuro((((150*(b - 1)) + 1)):end)};
        else
            % hippocampal
            hippTempDeltaStruct(b,1) = {hippDeltaNeuro((((150*(b - 1)) + 1)):(150*b))};
            hippTempThetaStruct(b,1) = {hippThetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            hippTempAlphaStruct(b,1) = {hippAlphaNeuro((((150*(b - 1)) + 1)):(150*b))};
            hippTempBetaStruct(b,1) = {hippBetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            hippTempGammaStruct(b,1) = {hippGammaNeuro((((150*(b - 1)) + 1)):(150*b))};
            hippTempMUAStruct(b,1) = {hippMUANeuro((((150*(b - 1)) + 1)):(150*b))};
            % cortical
            cortTempDeltaStruct(b,1) = {cortDeltaNeuro((((150*(b - 1)) + 1)):(150*b))};
            cortTempThetaStruct(b,1) = {cortThetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            cortTempAlphaStruct(b,1) = {cortAlphaNeuro((((150*(b - 1)) + 1)):(150*b))};
            cortTempBetaStruct(b,1) = {cortBetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            cortTempGammaStruct(b,1) = {cortGammaNeuro((((150*(b - 1)) + 1)):(150*b))};
            cortTempMUAStruct(b,1) = {cortMUANeuro((((150*(b - 1)) + 1)):(150*b))};
        end
    end
    % save hippocampal data under MergedData file
    MergedData.sleep.parameters.hippocampalNeural.deltaBandPower = hippTempDeltaStruct;
    MergedData.sleep.parameters.hippocampalNeural.thetaBandPower = hippTempThetaStruct;
    MergedData.sleep.parameters.hippocampalNeural.alphaBandPower = hippTempThetaStruct;
    MergedData.sleep.parameters.hippocampalNeural.betaBandPower = hippTempBetaStruct;
    MergedData.sleep.parameters.hippocampalNeural.gammaBandPower = hippTempGammaStruct;
    MergedData.sleep.parameters.hippocampalNeural.muaPower = hippTempMUAStruct;
    % save cortical data under MergedData file
    MergedData.sleep.parameters.corticalNeural.deltaBandPower = cortTempDeltaStruct;
    MergedData.sleep.parameters.corticalNeural.thetaBandPower = cortTempThetaStruct;
    MergedData.sleep.parameters.corticalNeural.alphaBandPower = cortTempAlphaStruct;
    MergedData.sleep.parameters.corticalNeural.betaBandPower = cortTempBetaStruct;
    MergedData.sleep.parameters.corticalNeural.gammaBandPower = cortTempGammaStruct;
    MergedData.sleep.parameters.corticalNeural.muaPower = cortTempMUAStruct;
    
    %% BLOCK PURPOSE: Create folder for the Neural spectrogram data of each electrode
    trialDuration_sec = 900;   % sec
    offset = 2.5;   % sec
    binWidth = 5;   % sec
    T = round(SpecData.corticalNeural.fiveSec.T,1);
    F = SpecData.corticalNeural.fiveSec.F;
    specCort = SpecData.corticalNeural.fiveSec.normS;
    specHip = SpecData.hippocampalNeural.fiveSec.normS;
    freqFloor = floor(F);
    % delta
    deltaLow = freqFloor == 1;
    deltaHigh = freqFloor == 4;
    deltaLowStart = find(deltaLow,1,'first');
    deltaLowEnd = find(deltaHigh,1,'last');
    deltaSpecHip = specHip(deltaLowStart:deltaLowEnd,:);
    deltaSpecCort = specCort(deltaLowStart:deltaLowEnd,:);
    meanDeltaSpecHip = mean(deltaSpecHip,1);
    meanDeltaSpecCort = mean(deltaSpecCort,1);
    % theta
    thetaLow = freqFloor == 4;
    thetaHigh = freqFloor == 10;
    thetaLowStart = find(thetaLow,1,'first');
    thetaLowEnd = find(thetaHigh,1,'last');
    thetaSpecHip = specHip(thetaLowStart:thetaLowEnd,:);
    thetaSpecCort = specCort(thetaLowStart:thetaLowEnd,:);
    meanThetaSpecHip = mean(thetaSpecHip,1);
    meanThetaSpecCort = mean(thetaSpecCort,1);
    % alpha
    alphaLow = freqFloor == 10;
    alphaHigh = freqFloor == 13;
    alphaLowStart = find(alphaLow,1,'first');
    alphaLowEnd = find(alphaHigh,1,'last');
    alphaSpecHip = specHip(alphaLowStart:alphaLowEnd,:);
    alphaSpecCort = specCort(alphaLowStart:alphaLowEnd,:);
    meanAlphaSpecHip = mean(alphaSpecHip,1);
    meanAlphaSpecCort = mean(alphaSpecCort,1);
    % beta
    betaLow = freqFloor == 13;
    betaHigh = freqFloor == 30;
    betaLowStart = find(betaLow,1,'first');
    betaLowEnd = find(betaHigh,1,'last');
    betaSpecHip = specHip(betaLowStart:betaLowEnd,:);
    betaSpecCort = specCort(betaLowStart:betaLowEnd,:);
    meanBetaSpecHip = mean(betaSpecHip,1);
    meanBetaSpecCort = mean(betaSpecCort,1);
    % gamma
    gammaLow = freqFloor == 30;
    gammaHigh = freqFloor == 99;
    gammaLowStart = find(gammaLow,1,'first');
    gammaLowEnd = find(gammaHigh,1,'last');
    gammaSpecHip = specHip(gammaLowStart:gammaLowEnd,:);
    gammaSpecCort = specCort(gammaLowStart:gammaLowEnd,:);
    meanGammaSpecHip = mean(gammaSpecHip,1);
    meanGammaSpecCort = mean(gammaSpecCort,1);
    % Divide the neural signals into five second bins and put them in a cell array
    hippTempDeltaSpecStruct = cell(180,1);
    hippTempThetaSpecStruct = cell(180,1);
    hippTempAlphaSpecStruct = cell(180,1);
    hippTempBetaSpecStruct = cell(180,1);
    hippTempGammaSpecStruct = cell(180,1);
    cortTempDeltaSpecStruct = cell(180,1);
    cortTempThetaSpecStruct = cell(180,1);
    cortTempAlphaSpecStruct = cell(180,1);
    cortTempBetaSpecStruct = cell(180,1);
    cortTempGammaSpecStruct = cell(180,1);
    % loop through all samples across the 15 minutes in 5 second bins (180 total)
    for c = 1:180
        if c == 1
            startTime = offset;
            startTime_index = find(T == startTime);
            endTime = 5;
            [~,endTime_index] = min(abs(T - endTime));
            % hippocampal
            hippTempDeltaSpecStruct{c,1} = {meanDeltaSpecHip(startTime_index:endTime_index)};
            hippTempThetaSpecStruct{c,1} = {meanThetaSpecHip(startTime_index:endTime_index)};
            hippTempAlphaSpecStruct{c,1} = {meanAlphaSpecHip(startTime_index:endTime_index)};
            hippTempBetaSpecStruct{c,1} = {meanBetaSpecHip(startTime_index:endTime_index)};
            hippTempGammaSpecStruct{c,1} = {meanGammaSpecHip(startTime_index:endTime_index)};
            % cortical
            cortTempDeltaSpecStruct{c,1} = {meanDeltaSpecCort(startTime_index:endTime_index)};
            cortTempThetaSpecStruct{c,1} = {meanThetaSpecCort(startTime_index:endTime_index)};
            cortTempAlphaSpecStruct{c,1} = {meanAlphaSpecCort(startTime_index:endTime_index)};
            cortTempBetaSpecStruct{c,1} = {meanBetaSpecCort(startTime_index:endTime_index)};
            cortTempGammaSpecStruct{c,1} = {meanGammaSpecCort(startTime_index:endTime_index)};
        elseif c == 180
            startTime = trialDuration_sec - 5;
            [~,startTime_index] = min(abs(T - startTime));
            endTime = trialDuration_sec - offset;
            [~,endTime_index] = min(abs(T - endTime));
            % hippocampal
            hippTempDeltaSpecStruct{c,1} = {meanDeltaSpecHip(startTime_index:endTime_index)};
            hippTempThetaSpecStruct{c,1} = {meanThetaSpecHip(startTime_index:endTime_index)};
            hippTempAlphaSpecStruct{c,1} = {meanAlphaSpecHip(startTime_index:endTime_index)};
            hippTempBetaSpecStruct{c,1} = {meanBetaSpecHip(startTime_index:endTime_index)};
            hippTempGammaSpecStruct{c,1} = {meanGammaSpecHip(startTime_index:endTime_index)};
            % cortical
            cortTempDeltaSpecStruct{c,1} = {meanDeltaSpecCort(startTime_index:endTime_index)};
            cortTempThetaSpecStruct{c,1} = {meanThetaSpecCort(startTime_index:endTime_index)};
            cortTempAlphaSpecStruct{c,1} = {meanAlphaSpecCort(startTime_index:endTime_index)};
            cortTempBetaSpecStruct{c,1} = {meanBetaSpecCort(startTime_index:endTime_index)};
            cortTempGammaSpecStruct{c,1} = {meanGammaSpecCort(startTime_index:endTime_index)};
        else
            startTime = binWidth*(c - 1);
            [~,startTime_index] = min(abs(T - startTime));
            endTime = binWidth*c;
            [~,endTime_index] = min(abs(T - endTime));
            % hippocampal
            hippTempDeltaSpecStruct{c,1} = {meanDeltaSpecHip(startTime_index + 1:endTime_index + 1)};
            hippTempThetaSpecStruct{c,1} = {meanThetaSpecHip(startTime_index + 1:endTime_index + 1)};
            hippTempAlphaSpecStruct{c,1} = {meanAlphaSpecHip(startTime_index + 1:endTime_index + 1)};
            hippTempBetaSpecStruct{c,1} = {meanBetaSpecHip(startTime_index + 1:endTime_index + 1)};
            hippTempGammaSpecStruct{c,1} = {meanGammaSpecHip(startTime_index + 1:endTime_index + 1)};
            % cortical
            cortTempDeltaSpecStruct{c,1} = {meanDeltaSpecCort(startTime_index + 1:endTime_index + 1)};
            cortTempThetaSpecStruct{c,1} = {meanThetaSpecCort(startTime_index + 1:endTime_index + 1)};
            cortTempAlphaSpecStruct{c,1} = {meanAlphaSpecCort(startTime_index + 1:endTime_index + 1)};
            cortTempBetaSpecStruct{c,1} = {meanBetaSpecCort(startTime_index + 1:endTime_index + 1)};
            cortTempGammaSpecStruct{c,1} = {meanGammaSpecCort(startTime_index + 1:endTime_index + 1)};
        end
    end
    % save hippocampal data under MergedData file
    MergedData.sleep.parameters.hippocampalNeural.specDeltaBandPower = hippTempDeltaSpecStruct;
    MergedData.sleep.parameters.hippocampalNeural.specThetaBandPower = hippTempThetaSpecStruct;
    MergedData.sleep.parameters.hippocampalNeural.specAlphaBandPower = hippTempAlphaSpecStruct;
    MergedData.sleep.parameters.hippocampalNeural.specBetaBandPower = hippTempBetaSpecStruct;
    MergedData.sleep.parameters.hippocampalNeural.specGammaBandPower = hippTempGammaSpecStruct;
    % save cortical data under MergedData file
    MergedData.sleep.parameters.corticalNeural.specDeltaBandPower = cortTempDeltaSpecStruct;
    MergedData.sleep.parameters.corticalNeural.specThetaBandPower = cortTempThetaSpecStruct;
    MergedData.sleep.parameters.corticalNeural.specAlphaBandPower = cortTempAlphaSpecStruct;
    MergedData.sleep.parameters.corticalNeural.specBetaBandPower = cortTempBetaSpecStruct;
    MergedData.sleep.parameters.corticalNeural.specGammaBandPower = cortTempGammaSpecStruct;
    
    %% BLOCK PURPOSE: Create folder for binarized whisking and binarized force sensor
    binWhiskerAngle = MergedData.data.binWhiskerAngle;
    binForceSensor = MergedData.data.binForceSensorL;
    whiskerAngle = MergedData.data.whiskerAngle;
    whiskerAcceleration = diff(whiskerAngle,2);
    % Find the number of whiskerBins due to frame drops.
    whiskerBinNumber = ceil(length(binWhiskerAngle)/150);
    % Divide the signal into five second bins and put them in a cell array
    tempWhiskerStruct = cell(whiskerBinNumber,1);
    tempWhiskerAccelStruct = cell(whiskerBinNumber,1);
    tempBinWhiskerStruct = cell(whiskerBinNumber,1);
    tempForceStruct = cell(whiskerBinNumber,1);
    for whiskerBins = 1:whiskerBinNumber
        if whiskerBins == 1
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle(whiskerBins:150)};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration(whiskerBins:150)};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle(whiskerBins:150)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor(whiskerBins:150)};
        elseif whiskerBins == whiskerBinNumber
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle((((150*(whiskerBins-1)) + 1)):end)};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration((((150*(whiskerBins-1)) + 1)):end)};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((150*(whiskerBins-1)) + 1)):end)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((150*(whiskerBins-1)) + 1)):end)};
        else
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
        end
    end
    % save whisker and force sensor data under MergedData file
    MergedData.sleep.parameters.whiskerAngle = tempWhiskerStruct;
    MergedData.sleep.parameters.whiskerAcceleration = tempWhiskerAccelStruct;
    MergedData.sleep.parameters.binWhiskerAngle = tempBinWhiskerStruct;
    MergedData.sleep.parameters.binForceSensor = tempForceStruct;
    
    %% Create folder for the EMG
    EMG = MergedData.data.EMG.data;
    normEMG = EMG - RestingBaselines.(baselineType).EMG.data.(strDay);  
    tempEMGStruct = cell(180,1);
    for EMGBins = 1:180
        if EMGBins == 1
            tempEMGStruct(EMGBins,1) = {normEMG(EMGBins:150)};
        else
            tempEMGStruct(EMGBins,1) = {normEMG((((150*(EMGBins-1)) + 1)):(150*EMGBins))};
        end
    end
    % save EMG data under MergedData file
    MergedData.sleep.parameters.EMG = tempEMGStruct;
    
    %% BLOCK PURPOSE: Create folder for the left and right CBV data
    vesselDiameter = MergedData.data.vesselDiameter.data;
    normVesselDiameter = (vesselDiameter - RestingBaselines.(baselineType).vesselDiameter.data.(vesselID).(strDay))/RestingBaselines.(baselineType).vesselDiameter.data.(vesselID).(strDay);
    tempVesselDiameterStruct = cell(180,1);
    for vesselBins = 1:180
        if vesselBins == 1
            tempVesselDiameterStruct(vesselBins,1) = {normVesselDiameter(vesselBins:25)};  % Samples 1 to 150
        else
            tempVesselDiameterStruct(vesselBins,1) = {normVesselDiameter((((25*(vesselBins-1)) + 1)):(25*vesselBins))};  % Samples 151 to 300, etc...
        end
    end
    % save hemodynamic data under MergedData file
    MergedData.sleep.parameters.vesselDiameter.data = tempVesselDiameterStruct;
    save(mergedDataFileID,'MergedData');
end

end
