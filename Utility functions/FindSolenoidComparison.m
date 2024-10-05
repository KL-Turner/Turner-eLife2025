function [comparison] = FindSolenoidComparison(hemisphere,solenoid)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
if strcmp(solenoid,'AudSol') == true
    comparison = 'aud';
elseif strcmp(solenoid,'LPadSol')
    if any(strcmp(hemisphere,{'LH','fLH'}))
        comparison = 'ipsi';
    else
        comparison = 'contra';
    end
elseif strcmp(solenoid,'RPadSol')
    if any(strcmp(hemisphere,{'RH','fRH'}))
        comparison = 'ipsi';
    else
        comparison = 'contra';
    end
end