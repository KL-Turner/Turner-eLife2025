%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
lineScanFs = MScanData.data.bloodFlow.Fs;
desiredFs = 5;
originalVelocity = abs(MScanData.data.bloodFlow.fixedVelocity);
% lowpass filter signal below Nyquist frequency of downsample rate;
[z,p,k] = butter(3,(0.5*floor(desiredFs))/(0.5*lineScanFs),'low');
[sos,g] = zp2sos(z,p,k);
filtOriginalVelocity = filtfilt(sos,g,originalVelocity);
% generate sine wave for resampling
nonIntegerSineWave = dsp.SineWave(1,desiredFs,1,'SampleRate',lineScanFs,'SamplesPerFrame',length(originalVelocity));
triggerWave = nonIntegerSineWave();
waveDiff = diff(triggerWave);
downSlope = waveDiff < 0;
edgeFind = diff(downSlope);
downsampleInds = edgeFind == 1;
downsampledVelocity = filtOriginalVelocity(:,downsampleInds);
% check figure
figure
p1 = plot((1:length(originalVelocity))/lineScanFs,originalVelocity,'r');
hold on;
p2 = plot((1:length(downsampledVelocity))/desiredFs,downsampledVelocity,'b');
xlabel('Time (s)')
ylabel('Velocity (\muM/sec)')
legend([p1,p2],'Original signal','Downsampled signal')