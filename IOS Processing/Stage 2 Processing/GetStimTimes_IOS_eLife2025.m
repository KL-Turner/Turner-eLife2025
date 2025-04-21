function [stimTimes] = GetStimTimes_IOS_eLife2025(ProcData)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% get solenoid times
stimNames = fieldnames(ProcData.data.stimulations);
stimList = cell(1,length(stimNames));
for sN = 1:length(stimNames)
    stimList{sN} = ProcData.data.stimulations.(stimNames{sN});
end
stimTimes = cell2mat(stimList);