function [r] = CalculateR_IOS_eLife2025(pred,act)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%
% Purpose: 
%________________________________________________________________________________________________________________________

% calculate correlation coefficient
coeffMat = corrcoef(pred,act);
r = coeffMat(2,1);

end
