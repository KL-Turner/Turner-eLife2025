function [procDataFileIDs] = CategorizeData_IOS()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% character list of all ProcData files in the directory from ProcessRawDataFiles_IOS.m
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    % load and Setup
    load(procDataFileID)
    if isfield(ProcData,'flags') == false
        disp(['Categorizing behavioral data for file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
        % condense pulsed stimulations
        stimulationFields = fieldnames(ProcData.data.stimulations);
        for vv = 1:length(stimulationFields)
            stimulationTimes = ProcData.data.stimulations.(stimulationFields{vv,1});
            condensedStimulationTimes = [];
            cc = 1;
            if isempty(stimulationTimes) == false
                for bb = 1:length(stimulationTimes)
                    if bb == 1
                        condensedStimulationTimes(1,bb) = stimulationTimes(1,bb);
                        cc = cc + 1;
                    else
                        timeDifference = stimulationTimes(1,bb) - stimulationTimes(1,bb - 1);
                        if timeDifference > 1 % remove stimulations that are closer than 1 second to the previous
                            condensedStimulationTimes(1,cc) = stimulationTimes(1,bb);
                            cc = cc + 1;
                        end
                    end
                end
                ProcData.data.stimulations.(stimulationFields{vv,1}) = condensedStimulationTimes;
            end
        end
        whiskerSamplingRate = ProcData.notes.dsFs;
        linkThresh = 0.5; % link events < 0.5 seconds apart
        breakThresh = 0;
        modBinWhiskers = ProcData.data.whiskerAngle.binarization;
        % link the binarized whisking
        binWhiskers = LinkBinaryEvents_IOS(gt(modBinWhiskers,0),[linkThresh breakThresh]*whiskerSamplingRate);
        % handle edge conditions
        if binWhiskers(1) == 0 && binWhiskers(2) == 1
            binWhiskers(1) = 1;
        elseif binWhiskers(1) == 1 && binWhiskers(2) == 0
            binWhiskers(1) = 0;
        end
        % handle edge conditions
        if binWhiskers(end) == 0 && binWhiskers(end - 1) == 1
            binWhiskers(end) = 1;
        elseif binWhiskers(end) == 1 && binWhiskers(end - 1) == 0
            binWhiskers(end) = 0;
        end
        % retrieve details on whisking events
        [ProcData.flags.whisk] = GetWhiskingData_IOS(ProcData,binWhiskers);
        % retrieve details on stiming events
        [ProcData.flags.stim] = GetStimData_IOS(ProcData);
        % identify and separate resting data
        [ProcData.flags.rest] = GetRestData_IOS(ProcData);
        % Save ProcData structure
        save(procDataFileID,'ProcData');
    else
        disp(['Categorization for file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ') already analyzed.']); disp(' ')
    end
end