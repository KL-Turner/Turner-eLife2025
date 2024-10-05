function [Results_GFP] = AnalyzeHemoGFPRelationship_EGFP_nNOS(animalID,group,set,rootFolder,delim,Results_GFP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
data.dummyCheck = 1;
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID)
    [~,fileDate,~] = GetFileInfo_IOS_nNOS(procDataFileID);
    strDay = ConvertDate_IOS_nNOS(fileDate);
    LH_blueChannel = ProcData.data.GCaMP7s.LH;
    RH_blueChannel = ProcData.data.GCaMP7s.RH;
    LH_greenChannel = ProcData.data.CBV.LH;
    RH_greenChannel = ProcData.data.CBV.RH;
    if isfield(data,'strDay') == false
        data.(strDay) = [];
        data.(strDay).LH_blue = LH_blueChannel;
        data.(strDay).RH_blue = RH_blueChannel;
        data.(strDay).LH_green = LH_greenChannel;
        data.(strDay).RH_green = RH_greenChannel;
    else
        data.(strDay).LH_blue = cat(2,data.(strDay).LH_blueLH_blueChannel);
        data.(strDay).RH_blue = cat(2,data.(strDay).RH_blue,RH_blueChannel);
        data.(strDay).LH_green = cat(2,data.(strDay).LH_green,LH_greenChannel);
        data.(strDay).RH_green = cat(2,data.(strDay).RH_green,RH_greenChannel);
    end
end
uniqueDays = fieldnames(data);
finalData.LH_blue = []; finalData.RH_blue = []; finalData.LH_green = []; finalData.RH_green = [];
for aa = 1:length(uniqueDays)
    finalData.LH_blue = cat(2,finalData.LH_blue,(data.(strDay).LH_blue - mean(data.(strDay).LH_blue))./mean(data.(strDay).LH_blue));
    finalData.RH_blue =  cat(2,finalData.RH_blue,(data.(strDay).RH_blue - mean(data.(strDay).RH_blue))./mean(data.(strDay).RH_blue));
    finalData.LH_green =  cat(2,finalData.LH_green,(data.(strDay).LH_green - mean(data.(strDay).LH_green))./mean(data.(strDay).LH_green));
    finalData.RH_green =  cat(2,finalData.RH_green,(data.(strDay).RH_green - mean(data.(strDay).RH_green))./mean(data.(strDay).RH_green));
end
Results_GFP.(group).(animalID).LH.blue = finalData.LH_blue;
Results_GFP.(group).(animalID).RH.blue = finalData.RH_blue;
Results_GFP.(group).(animalID).LH.green = finalData.LH_green;
Results_GFP.(group).(animalID).RH.green = finalData.RH_green;
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_GFP.mat','Results_GFP')
cd([rootFolder delim 'Data'])