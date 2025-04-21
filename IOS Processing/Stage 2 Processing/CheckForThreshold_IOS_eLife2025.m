function [ok] = CheckForThreshold_IOS_eLife2025(sfield,animal)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% begin Check
ok = 0;
if exist([animal '_Thresholds.mat'],'file') == 2
    load([animal '_Thresholds.mat']);
    if isfield(Thresholds,sfield)
        ok = 1;
    end
end