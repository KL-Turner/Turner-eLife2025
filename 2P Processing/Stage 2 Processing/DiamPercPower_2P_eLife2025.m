function [diamPerc,S,f] = DiamPercPower_2P_eLife2025(rawDiameter,baseDiameter,samplingRate)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: Calculates the percent of diameter power.
%________________________________________________________________________________________________________________________

L_Diam = length(rawDiameter);
duration = L_Diam/samplingRate;
the_fs = 1/samplingRate;
% change raw diameter to diameter percentage with a median filter and a low pass filter
d_Amp = (rawDiameter/baseDiameter - 1)*100;
diamMedf = medfilt1(d_Amp,5);   % 5-point median filter
[B,A] = butter(3,3/samplingRate,'low');   % change from 6 to 3
diamPerc = filtfilt(B,A,diamMedf);
% power spectrum multitaper
NW = floor(max(1,the_fs*duration/2));
params.Fs = samplingRate;
params.tapers = [NW,2*NW - 1];
params.err = [1,.01];
params.fpass = [0.05,3];
[S,f,~] = mtspectrumc(diamPerc - mean(diamPerc),params);

end
