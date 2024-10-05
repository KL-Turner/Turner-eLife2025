function [r2] = CalculateRsquared_HRF2020(pred,act)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%
% Purpose: 
%________________________________________________________________________________________________________________________

% % error Variance
% SSE = sum((act - pred).^2);
% % total Variance
% SST = sum((act - (ones(size(act,1),1)*mean(act))).^2);
% % check that the sum of the residuals is small compared to SSE + SSR
% r2 = ones(size(SSE)) - SSE./SST;


mdl = fitlm(act,pred);
r2 = mdl.Rsquared.Ordinary;

end
