function [MScanData] = GetVelocityRadon_2P_eLife2025(fileID,angleSpan,angleResolution,maxVelocity,interleaved,MScanData)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

frames = [MScanData.notes.startFrame,MScanData.notes.endFrame];
theX = [MScanData.notes.xStart,MScanData.notes.xStop];
windowSizeA = round(MScanData.notes.lineRate*MScanData.notes.frameTime);
% adapt windowsize to interleaved/not interleaved
if interleaved == 1
    windowSize = 8*round(windowSizeA/8);
else
    windowSize = 4*round(windowSizeA/4);
end
MScanData.data.bloodFlow.windowsize = windowSize;
nFrames = frames(2) - frames(1) + 1;
stepSize = .25*windowSize;
nLines = nFrames*(MScanData.notes.ySize);
nSteps = floor(nLines/stepSize) - 3;
angles = (1:180);
anglesBaseline = 1:15:180;
anglesAdaptive = -angleSpan:1:angleSpan; % angle range to look over
anglesFine = -1.75:angleResolution:1.75; % fine grained search
adaptiveAngles = length(anglesAdaptive);
thetaLow = 5; % boundaries for for adaptive theta
thetaHigh = adaptiveAngles - 5;
spreadMatrix = zeros(1,length(angles));
spreadMatrixAdaptive = zeros(2,-nSteps,length(anglesAdaptive));
baselineVar = zeros(2,nSteps,length(anglesBaseline));
spreadMatrixFine = zeros(2,nSteps,length(anglesFine));
thetas = zeros(2,nSteps); % the angle of the lines
% data caching
dataTemp = LoadTiffConcatenate_2P_eLife2025([fileID '.TIF'],frames);
dataCache = dataTemp(:,theX(1):MScanData.notes.decimate:theX(2) - 1);
MScanData.data.bloodFlow.mean_BG = (mean(dataCache));
useLines = 1:windowSize;
if interleaved == 1
    dataHold = zeros(2,windowSize/2,length(theX(1):MScanData.notes.decimate:theX(2) - 1));
    dataHold(1,:,:) = double(dataCache(1:2:windowSize,1:MScanData.notes.decimate:end));
    dataHold(2,:,:) = double(dataCache(2:2:windowSize,end:-MScanData.notes.decimate:1));
else
    dataHold = zeros(1,windowSize,length(theX(1):MScanData.notes.decimate:theX(2) - 1));
    dataHold(1,:,:) = double(dataCache(1:windowSize,1:MScanData.notes.decimate:end));
end
lineCounter = 4*stepSize; % this is the end of the window
dataCacheMeanForward = mean(dataCache);
if interleaved == 1
    dataCacheMeanBackward = (fliplr(dataCacheMeanForward));
end
if interleaved == 1
    for kk = 1:nSteps
        frameLineMarker = lineCounter; % find the new place in the frame
        dataHold(1,:,:) = double(dataCache(frameLineMarker - 4*stepSize + 1:2:frameLineMarker - 1,1:MScanData.notes.decimate:end));
        dataHold(2,:,:) = double(dataCache(frameLineMarker + 2 - 4*stepSize:2:frameLineMarker,end:-MScanData.notes.decimate:1));
        useLines=useLines+stepSize;
        % take out the local mean of the window
        theT(kk) = 1 + (kk - 1)*stepSize + windowSize/2;
        npointsDecimated = size(dataHold,3);
        for n = 1:npointsDecimated
            dataHold_ms(1,:,n) = dataHold(1,:,n) - dataCacheMeanForward(n);
            dataHold_ms(2,:,n) = dataHold(2,:,n) - dataCacheMeanBackward(n);
        end
        dataHold_ms = dataHold_ms - mean(dataHold_ms(:));
        if kk == 1
            radonHold(1,:,:) = (radon(squeeze(dataHold_ms(1,:,:)),angles));
            radonHold(2,:,:) = radon(squeeze(dataHold_ms(2,:,:)),angles);
            spreadMatrix(1,:) = var(squeeze(radonHold(1,:,:)));
            spreadMatrix(2,:) = var(squeeze(radonHold(2,:,:)));
            [~,theTheta(1,:)] = max(spreadMatrix(1,kk));
            [~,theTheta(2,:)] = max(spreadMatrix(2,kk));
            thetas(1,kk) = angles(theTheta(1,:));
            thetas(2,kk) = angles(theTheta(2,:));
            radonHoldAdaptive(1,:,:) = radon(squeeze(dataHold_ms(1,:,:)),-thetas(1,kk) + anglesAdaptive);
            radonHoldAdaptive(2,:,:) = radon(squeeze(dataHold_ms(2,:,:)),-thetas(2,kk) + anglesAdaptive);
            spreadMatrixAdaptive(1,kk,1) = max(var(squeeze(radonHoldAdaptive(1,:,:))));
            spreadMatrixAdaptive(2,kk,1) = max(var(squeeze(radonHoldAdaptive(2,:,:))));
            radonHoldFine(1,:,:) = radon(squeeze(dataHold_ms(1,:,:) - mean(dataHold_ms(:))),thetas(1,kk) + anglesFine);
            radonHoldFine(2,:,:) = radon(squeeze(dataHold_ms(2,:,:) - mean(dataHold_ms(:))),thetas(2,kk) + anglesFine);
            spreadMatrixFine(1,kk,:) = squeeze(var(squeeze(radonHoldFine(1,:,:))));
            spreadMatrixFine(2,kk,:) = squeeze(var(squeeze(radonHoldFine(2,:,:))));
            [~,theTheta(1,:)] = max(squeeze(spreadMatrixFine(1,kk,:)));
            [~,theTheta(2,:)] = max(squeeze(spreadMatrixFine(2,kk,:)));
            thetas(1,kk) = thetas(1,kk) + anglesFine(theTheta(1,:));
            thetas(2,kk) = thetas(2,kk) + anglesFine(theTheta(2,:));
            baselineVar(1,kk,:) = var(radon(squeeze(dataHold_ms(1,:,:)),anglesBaseline));
            baselineVar(2,kk,:) = var(radon(squeeze(dataHold_ms(2,:,:)),anglesBaseline));
        else
            thetas(1,kk) = thetas(1,kk - 1);
            thetas(2,kk) = thetas(2,kk - 1);
            radonHoldAdaptive(1,:,:) = radon(squeeze(dataHold_ms(1,:,:)),thetas(1,kk - 1) + anglesAdaptive);
            radonHoldAdaptive(2,:,:) = radon(squeeze(dataHold_ms(2,:,:)),thetas(2,kk - 1) + anglesAdaptive);
            spreadMatrixAdaptive(1,kk,:) = var(squeeze(radonHoldAdaptive(1,:,:)));
            spreadMatrixAdaptive(2,kk,:) = var(squeeze(radonHoldAdaptive(2,:,:)));
            [~,theTheta(1,:)] = max(spreadMatrixAdaptive(1,kk,:));
            [~,theTheta(2,:)] = max(spreadMatrixAdaptive(2,kk,:));
            % figure out the average radon of the section for calculating S/N
            baselineVar(1,kk,:) = var(radon(squeeze(dataHold_ms(1,:,:)),anglesBaseline));
            baselineVar(2,kk,:) = var(radon(squeeze(dataHold_ms(2,:,:)),anglesBaseline));
            % if the peak of the variance is at an extreme, redo for angles centered
            while ((theTheta(1,:)<= thetaLow)||(theTheta(1,:) >= thetaHigh))
                thetas(1,kk) = thetas(1,kk) + anglesAdaptive(theTheta(1,:));
                radonHoldAdaptive(1,:,:) = radon(squeeze(dataHold_ms(1,:,:)),thetas(1,kk) + anglesAdaptive);
                spreadMatrixAdaptive(1,kk,:) = var(squeeze(radonHoldAdaptive(1,:,:)));
                [~,theTheta(1,:)] = max(squeeze(spreadMatrixAdaptive(1,kk,:)));
            end
            while ((theTheta(2,:) <= thetaLow)||(theTheta(2,:) >= thetaHigh))
                thetas(2,kk) = thetas(2,kk) + anglesAdaptive(theTheta(2,:));
                radonHoldAdaptive(2,:,:) = radon(squeeze(dataHold_ms(2,:,:)),thetas(2,kk) + anglesAdaptive);
                spreadMatrixAdaptive(2,kk,:) = var(squeeze(radonHoldAdaptive(2,:,:)));
                [~,theTheta(2,:)] = max(squeeze(spreadMatrixAdaptive(2,kk,:)));
            end
            baselineVar(1,kk,:) = var(radon(squeeze(dataHold_ms(1,:,:)),anglesBaseline));
            baselineVar(2,kk,:) = var(radon(squeeze(dataHold_ms(2,:,:)),anglesBaseline));
            thetas(1,kk) = thetas(1,kk) + anglesAdaptive(theTheta(1,:));
            radonHoldFine(1,:,:) = radon(squeeze(dataHold_ms(1,:,:)),thetas(1,kk) + anglesFine);
            spreadMatrixFine(1,kk,:) = var(squeeze(radonHoldFine(1,:,:)));
            [~,theTheta(1,:)] = max(spreadMatrixFine(1,kk,:));
            thetas(1,kk) = thetas(1,kk) + anglesFine(theTheta(1,:));
            thetas(2,kk) = thetas(2,kk) + anglesAdaptive(theTheta(2,:));
            radonHoldFine(2,:,:) = radon(squeeze(dataHold_ms(2,:,:)),thetas(2,kk) + anglesFine);
            spreadMatrixFine(2,kk,:)=var(squeeze(radonHoldFine(2,:,:)));
            [~,theTheta(2,:)] = max(spreadMatrixFine(2,kk,:));
            thetas(2,kk) = thetas(2,kk) + anglesFine(theTheta(2,:));
        end
        lineCounter = lineCounter + stepSize;
    end
    MScanData.data.bloodFlow.thetas = thetas;
    MScanData.data.bloodFlow.unscaledVelocity = (cotd(thetas - 90));
    MScanData.data.bloodFlow.spreadMatrixAdaptive = spreadMatrixAdaptive;
    MScanData.data.bloodFlow.baselineVar = baselineVar;
    % need to divide by 2 because of splitting due to inerleaving
    if interleaved == 1
        MScanData.data.bloodFlow.v(1,:) = MScanData.notes.decimate*MScanData.notes.tFactor*MScanData.notes.xFactor*(cotd(thetas(1,:) - 90))/2;
        MScanData.data.bloodFlow.v(2,:) = MScanData.notes.decimate*MScanData.notes.tFactor*MScanData.notes.xFactor*(cotd(thetas(2,:) - 90))/2;
    else
        MScanData.data.bloodFlow.v(1,:) = MScanData.notes.decimate*MScanData.notes.tFactor*MScanData.notes.xFactor*(cotd(thetas(1,:) - 90));
        MScanData.data.bloodFlow.v(2,:) = MScanData.notes.decimate*MScanData.notes.tFactor*MScanData.notes.xFactor*(cotd(thetas(1,:) - 90));
    end
    MScanData.data.bloodFlow.velocity2(1,:) = VelocityCleanUp_2P_eLife2025(MScanData.data.bloodFlow.v(1,:),maxVelocity);
    if interleaved == 1
        MScanData.data.bloodFlow.velocity2(2,:) = VelocityCleanUp_2P_eLife2025(MScanData.data.bloodFlow.v(2,:),maxVelocity);
    else
        MScanData.data.bloodFlow.velocity2(2,:) = -(MScanData.data.bloodFlow.velocity2(1,:));
    end
    vHold = MScanData.data.bloodFlow.velocity2;
    vOut = 0.5*(vHold(1,:) - vHold(2,:)); % take the average of the two directions
    MScanData.data.bloodFlow.vOut = vOut;
    MScanData.data.bloodFlow.theT = theT/(MScanData.notes.lineRate);
    MScanData.data.bloodFlow.Fs = 4/MScanData.notes.frameTime;
    [MScanData.data.bloodFlow.fixedVelocity] = VelocityCleanUpSTDOutliers_2P_eLife2025(MScanData.data.bloodFlow.vOut,3); % clear outliers > 3 standard deviations
    try
        MScanData.data.bloodFlow.Image = dataCache;
    catch
        MScanData.data.bloodFlow.Image = zeros(MScanData.notes.num_frames*(MScanData.notes.Frame_Height),1);
    end
    MScanData.data.bloodFlow.separability = max(squeeze(mean(MScanData.data.bloodFlow.spreadMatrixAdaptive,1))')./mean(squeeze(mean(MScanData.data.bloodFlow.baselineVar,1))');%#ok<UDIM> %
    % figures
    lineScanVel = figure;
    sgtitle([MScanData.notes.animalID ' ' MScanData.notes.date ' ' MScanData.notes.imageID])
    subplot(3,1,1)
    p1 = plot(MScanData.data.bloodFlow.theT,MScanData.data.bloodFlow.velocity2(1,:));
    hold on
    p2 = plot(MScanData.data.bloodFlow.theT,MScanData.data.bloodFlow.velocity2(2,:));
    p3 = plot(MScanData.data.bloodFlow.theT,MScanData.data.bloodFlow.vOut);
    title('Velocities');
    xlabel('Time (s)');
    ylabel('Velocity (um/s)');
    legend([p1,p2,p3],'Forward','Backward','vOut');
    set(gca,'box','off')
    subplot(4,1,2)
    plot(MScanData.data.bloodFlow.theT,MScanData.data.bloodFlow.fixedVelocity);
    title('Fixed Velocity')
    xlabel('Time (s)')
    ylabel('Velocity (um/s)')
    set(gca,'box','off')
    subplot(4,1,3)
    plot(MScanData.data.bloodFlow.theT,MScanData.data.bloodFlow.thetas(1,:));
    hold on
    plot(MScanData.data.bloodFlow.theT,MScanData.data.bloodFlow.thetas(2,:));
    title('Thetas')
    xlabel('Time (s)')
    ylabel('Theta (degrees)')
    legend('Theta1','Theta2')
    set(gca,'box','off')
    subplot(4,1,4)
    plot(MScanData.data.bloodFlow.theT,MScanData.data.bloodFlow.separability);
    title('Separability (> 3 = good)')
    xlabel('Time (s)')
    ylabel('Separability')
    set(gca,'box','off')
    % save path
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Thetas and Separability/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(lineScanVel,[dirpath MScanData.notes.animalID '_' MScanData.notes.date '_' MScanData.notes.imageID '_VelocityThetaSeparability']);
    close(lineScanVel)
end

end
