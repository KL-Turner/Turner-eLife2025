function out = genFigure_individual_SAP(fileList)
HR.LTA = [];

numLTA = 0;

HbT.PC.LTA = [];

HbD.PC.LTA = [];

% resting state parameters

RestVar.PC.HbT = [];

RestVar.PC.HbD = [];

RestVar.PC.HR = [];


for file_idx = 1:numel(fileList)
    filepath = fileList{file_idx};
    data = load(filepath);
    parts = strsplit(filepath, filesep);
    animal = parts{end-1};
    file_name = parts{end};
    
    
    HR.LTA = [HR.LTA;data.LTA_HR;]; % Heart rate
    
    HbT.PC.LTA = [HbT.PC.LTA; data.LTA_PC_HbO.RH + data.LTA_PC_HbR.RH;];
    HbD.PC.LTA = [HbD.PC.LTA; data.LTA_PC_HbO.RH - data.LTA_PC_HbR.RH;];
    numLTA = numLTA + data.LTA_STDist_num;
    
    RestVar.PC.HR = [RestVar.PC.HR; data.RestVarInd.HR;];
    RestVar.PC.HbT = [RestVar.PC.HbT; data.RestVarInd.HbT.PC;];
    RestVar.PC.HbD = [RestVar.PC.HbD; data.RestVarInd.HbD.PC;];
    
end

%% Output data
out.HR = HR;


out.numLTA = numLTA;

out.HbT = HbT;
out.HbD = HbD;

out.RestVar = RestVar;



    
  
end
