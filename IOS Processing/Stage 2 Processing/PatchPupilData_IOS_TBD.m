function [ProcData] = PatchPupilData_IOS(RawData,ProcData)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% expected number of frames based on trial duration and sampling rate
expectedSamples = RawData.notes.trialDuration_sec*RawData.notes.pupilCamSamplingRate;
droppedFrameIndex = str2double(RawData.notes.droppedPupilCamFrameIndex);
sampleDiff = expectedSamples - length(RawData.data.pupil.pupilMajor);
framesPerIndex = ceil(sampleDiff/length(droppedFrameIndex));
% patch NaN values
pupilArea = RawData.data.pupil.pupilArea;
pupilMinorAxis = RawData.data.pupil.pupilMinor;
pupilMajorAxis = RawData.data.pupil.pupilMajor;
pupilCentroidX = RawData.data.pupil.pupilCentroid(:,1)';
pupilCentroidY = RawData.data.pupil.pupilCentroid(:,2)';
% patch NaNs due to blinking
blinkNaNs = isnan(pupilArea);
[linkedBlinkIndex] = LinkBinaryEvents_IOS(gt(blinkNaNs,0),[RawData.notes.pupilCamSamplingRate,0]); % link greater than 1 second
% identify edges for interpolation
xx = 1;
edgeFoundA = false;
startEdgeA = [];
endEdgeA = [];
for aa = 1:length(linkedBlinkIndex)
    if edgeFoundA == false
        if linkedBlinkIndex(1,aa) == 1 && (aa < length(linkedBlinkIndex)) == true
            startEdgeA(xx,1) = aa;
            edgeFoundA = true;
        end
    elseif edgeFoundA == true
        if linkedBlinkIndex(1,aa) == 0
            endEdgeA(xx,1) = aa;
            edgeFoundA = false;
            xx = xx + 1;
        elseif (length(linkedBlinkIndex) == aa) == true && (linkedBlinkIndex(1,aa) == 1) == true
            endEdgeA(xx,1) = aa;
        end
    end
end
% fill from start:ending edges of rapid pupil fluctuations that weren't NaN
for aa = 1:length(startEdgeA)
    try
        pupilArea(startEdgeA(aa,1) - 2:endEdgeA(aa,1) + 2) = NaN;
        pupilMinorAxis(startEdgeA(aa,1) - 2:endEdgeA(aa,1) + 2) = NaN;
        pupilMajorAxis(startEdgeA(aa,1) - 2:endEdgeA(aa,1) + 2) = NaN;
        pupilCentroidX(startEdgeA(aa,1) - 2:endEdgeA(aa,1) + 2) = NaN;
        pupilCentroidY(startEdgeA(aa,1) - 2:endEdgeA(aa,1) + 2) = NaN;
        patchLength(aa,1) = (endEdgeA(aa,1) + 2) - (startEdgeA(aa,1) - 2);
    catch
        pupilArea(startEdgeA(aa,1):endEdgeA(aa,1)) = NaN;
        pupilMinorAxis(startEdgeA(aa,1):endEdgeA(aa,1)) = NaN;
        pupilMajorAxis(startEdgeA(aa,1):endEdgeA(aa,1)) = NaN;
        pupilCentroidX(startEdgeA(aa,1):endEdgeA(aa,1)) = NaN;
        pupilCentroidY(startEdgeA(aa,1):endEdgeA(aa,1)) = NaN;
        patchLength(aa,1) = endEdgeA(aa,1) - startEdgeA(aa,1);
    end
end
% patch NaN values with moving median filter
try
    patchedPupilArea = fillmissing(pupilArea,'movmedian',max(patchLength)*2);
    patchedMinorAxis = fillmissing(minorAxis,'movmedian',max(patchLength)*2);
    patchedMajorAxis = fillmissing(majorAxis,'movmedian',max(patchLength)*2);
    patchedCentroidX = fillmissing(pupilCentroidX,'movmedian',max(patchLength)*2);
    patchedCentroidY = fillmissing(pupilCentroidY,'movmedian',max(patchLength)*2);
catch
    patchedPupilArea = pupilArea;
    patchedMinorAxis = pupilMinorAxis;
    patchedMajorAxis = pupilMajorAxis;
    patchedCentroidX = pupilCentroidX;
    patchedCentroidY = pupilCentroidY;
end
%% patch sudden spikes
diffArea = abs(diff(patchedPupilArea));
% threshold for interpolation
threshold = 250;
diffIndex = diffArea > threshold;
[linkedDiffIndex] = LinkBinaryEvents_IOS(gt(diffIndex,0),[RawData.notes.pupilCamSamplingRate*2,0]);
% identify edges for interpolation
edgeFoundB = false;
xx = 1;
startEdgeB = [];
endEdgeB = [];
for aa = 1:length(linkedDiffIndex)
    if edgeFoundB == false
        if (linkedDiffIndex(1,aa) == 1) == true && (aa < length(linkedDiffIndex)) == true
            startEdgeB(xx,1) = aa;
            edgeFoundB = true;
        end
    elseif edgeFoundB == true
        if linkedDiffIndex(1,aa) == 0
            endEdgeB(xx,1) = aa;
            edgeFoundB = false;
            xx = xx + 1;
        elseif (length(linkedDiffIndex) == aa) == true && (linkedDiffIndex(1,aa) == 1) == true && edgeFoundB == true
            endEdgeB(xx,1) = aa;
        end
    end
end
% fill from start:ending edges of rapid pupil fluctuations that weren't NaN
for aa = 1:length(startEdgeB)
    try
        patchedPupilArea(startEdgeB(aa,1) - 2:endEdgeB(aa,1) + 2) = NaN;
        patchedMinorAxis(startEdgeB(aa,1) - 2:endEdgeB(aa,1) + 2) = NaN;
        patchedMajorAxis(startEdgeB(aa,1) - 2:endEdgeB(aa,1) + 2) = NaN;
        patchedCentroidX(startEdgeB(aa,1) - 2:endEdgeB(aa,1) + 2) = NaN;
        patchedCentroidY(startEdgeB(aa,1) - 2:endEdgeB(aa,1) + 2) = NaN;
        patchLength = (endEdgeB(aa,1) + 2) - (startEdgeB(aa,1) - 2);
    catch
        patchedPupilArea(startEdgeB(aa,1):endEdgeB(aa,1)) = NaN;
        patchedMinorAxis(startEdgeB(aa,1):endEdgeB(aa,1)) = NaN;
        patchedMajorAxis(startEdgeB(aa,1):endEdgeB(aa,1)) = NaN;
        patchedCentroidX(startEdgeB(aa,1):endEdgeB(aa,1)) = NaN;
        patchedCentroidY(startEdgeB(aa,1):endEdgeB(aa,1)) = NaN;
        patchLength = endEdgeB(aa,1) - startEdgeB(aa,1);
    end
    patchedPupilArea = fillmissing(patchedPupilArea,'movmedian',patchLength*2);
    patchedMinorAxis = fillmissing(patchedMinorAxis,'movmedian',patchLength*2);
    patchedMajorAxis = fillmissing(patchedMajorAxis,'movmedian',patchLength*2);
    patchedCentroidX = fillmissing(patchedCentroidX,'movmedian',patchLength*2);
    patchedCentroidY = fillmissing(patchedCentroidY,'movmedian',patchLength*2);
end
pupilArea = patchedPupilArea;
pupilMinorAxis = patchedMinorAxis;
pupilMajorAxis = patchedMajorAxis;
pupilCentroidX = patchedCentroidX;
pupilCentroidY = patchedCentroidY;
nanLogical = isnan(pupilArea);
nanIndex = find(nanLogical == 1);
if sum(nanIndex) > 1 && sum(nanIndex) < 1000
    while sum(nanIndex) >= 1
        pupilArea = fillmissing(pupilArea,'movmedian',3);
        pupilMinorAxis = fillmissing(pupilMinorAxis,'movmedian',3);
        pupilMajorAxis = fillmissing(pupilMajorAxis,'movmedian',3);
        pupilCentroidX = fillmissing(pupilCentroidX,'movmedian',3);
        pupilCentroidY = fillmissing(pupilCentroidY,'movmedian',3);
        nanIndex = isnan(pupilMajorAxis);
    end
end
%% patch missing frames now that NaN are gone
if ~isempty(droppedFrameIndex)
    % each dropped index
    for cc = 1:length(droppedFrameIndex)
        % for the first event, it's okay to start at the actual index
        if cc == 1
            leftEdge = droppedFrameIndex(1,cc);
        else
            % for all other dropped frames after the first, we need to correct for the fact that index is shifted right.
            leftEdge = droppedFrameIndex(1,cc) + ((cc - 1)*framesPerIndex);
        end
        % set the edges for the interpolation points. we want n number of samples between the two points,vthe left and
        % right edge values. This equates to having a 1/(dropped frames + 1) step size between the edges.
        rightEdge = leftEdge + 1;
        patchFrameInds = leftEdge:(1/(framesPerIndex + 1)):rightEdge;
        % concatenate the original data for the first index, then the new patched data for all subsequent
        % indeces. Take the values from 1:left edge, add in the new frames, then right edge to end.
        if cc == 1
            patchFrameVals_area = interp1(1:length(pupilArea),pupilArea,patchFrameInds); % linear interp
            patchFrameVals_minorAxis = interp1(1:length(pupilMinorAxis),pupilMinorAxis,patchFrameInds); % linear interp
            patchFrameVals_majorAxis = interp1(1:length(pupilMajorAxis),pupilMajorAxis,patchFrameInds); % linear interp
            patchFrameVals_centroidX = interp1(1:length(pupilCentroidX),pupilCentroidX,patchFrameInds); % linear interp
            patchFrameVals_centroidY = interp1(1:length(pupilCentroidY),pupilCentroidY,patchFrameInds); % linear interp
            snipPatchFrameVals_area = patchFrameVals_area(2:end - 1);
            snipPatchFrameVals_minorAxis = patchFrameVals_minorAxis(2:end - 1);
            snipPatchFrameVals_majorAxis = patchFrameVals_majorAxis(2:end - 1);
            snipPatchFrameVals_centroidX = patchFrameVals_centroidX(2:end - 1);
            snipPatchFrameVals_centroidY = patchFrameVals_centroidY(2:end - 1);
            try
                patchedArea = horzcat(pupilArea(1:leftEdge),snipPatchFrameVals_area,pupilArea(rightEdge:end));
                patchedMinorAxis = horzcat(pupilMinorAxis(1:leftEdge),snipPatchFrameVals_minorAxis,pupilMinorAxis(rightEdge:end));
                patchedMajorAxis = horzcat(pupilMajorAxis(1:leftEdge),snipPatchFrameVals_majorAxis,pupilMajorAxis(rightEdge:end));
                patchedCentroidX = horzcat(pupilCentroidX(1:leftEdge),snipPatchFrameVals_centroidX,pupilCentroidX(rightEdge:end));
                patchedCentroidY = horzcat(pupilCentroidY(1:leftEdge),snipPatchFrameVals_centroidY,pupilCentroidY(rightEdge:end));
            catch
                patchedArea = horzcat(pupilArea(1:end),snipPatchFrameVals_area);
                patchedMinorAxis = horzcat(pupilMinorAxis(1:end),snipPatchFrameVals_minorAxis);
                patchedMajorAxis = horzcat(pupilMajorAxis(1:end),snipPatchFrameVals_majorAxis);
                patchedCentroidX = horzcat(pupilCentroidX(1:end),snipPatchFrameVals_centroidX);
                patchedCentroidY = horzcat(pupilCentroidY(1:end),snipPatchFrameVals_centroidY);
            end
        else
            patchFrameVals_area = interp1(1:length(patchedArea),patchedArea,patchFrameInds); % linear interp
            patchFrameVals_minorAxis = interp1(1:length(patchedMinorAxis),patchedMinorAxis,patchFrameInds); % linear interp
            patchFrameVals_majorAxis = interp1(1:length(patchedMajorAxis),patchedMajorAxis,patchFrameInds); % linear interp
            patchFrameVals_centroidX = interp1(1:length(patchedCentroidX),patchedCentroidX,patchFrameInds); % linear interp
            patchFrameVals_centroidY = interp1(1:length(patchedCentroidY),patchedCentroidY,patchFrameInds); % linear interp
            snipPatchFrameVals_area = patchFrameVals_area(2:end - 1);
            snipPatchFrameVals_minorAxis = patchFrameVals_minorAxis(2:end - 1);
            snipPatchFrameVals_majorAxis = patchFrameVals_majorAxis(2:end - 1);
            snipPatchFrameVals_centroidX = patchFrameVals_centroidX(2:end - 1);
            snipPatchFrameVals_centroidY = patchFrameVals_centroidY(2:end - 1);
            patchedArea = horzcat(patchedArea(1:leftEdge),snipPatchFrameVals_area,patchedArea(rightEdge:end));
            patchedMinorAxis = horzcat(patchedMinorAxis(1:leftEdge),snipPatchFrameVals_minorAxis,patchedMinorAxis(rightEdge:end));
            patchedMajorAxis = horzcat(patchedMajorAxis(1:leftEdge),snipPatchFrameVals_majorAxis,patchedMajorAxis(rightEdge:end));
            patchedCentroidX = horzcat(patchedCentroidX(1:leftEdge),snipPatchFrameVals_centroidX,patchedCentroidX(rightEdge:end));
            patchedCentroidY = horzcat(patchedCentroidY(1:leftEdge),snipPatchFrameVals_centroidY,patchedCentroidY(rightEdge:end));
        end
    end
    % due to rounding up on the number of dropped frames per index, we have a few extra frames. Snip them off.
    patchedArea = patchedArea(1:expectedSamples);
    patchedMinorAxis = patchedMinorAxis(1:expectedSamples);
    patchedMajorAxis = patchedMajorAxis(1:expectedSamples);
    patchedCentroidX = patchedCentroidX(1:expectedSamples);
    patchedCentroidY = patchedCentroidY(1:expectedSamples);
else
    patchedArea = pupilArea(1:expectedSamples);
    patchedMinorAxis = pupilMinorAxis(1:expectedSamples);
    patchedMajorAxis = pupilMajorAxis(1:expectedSamples);
    patchedCentroidX = pupilCentroidX(1:expectedSamples);
    patchedCentroidY = pupilCentroidY(1:expectedSamples);
end
ProcData.data.pupil.area = patchedArea;
ProcData.data.pupil.minorAxis = patchedMinorAxis;
ProcData.data.pupil.majorAxis = patchedMajorAxis;
ProcData.data.pupil.centroidX = patchedCentroidX;
ProcData.data.pupil.centroidY = patchedCentroidY;
ProcData.data.blinks.index = RawData.data.pupil.blinkInds;