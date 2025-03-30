function [theCode]=GetFunctionCode_2P_nNOS(filename)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: this function opens the file and returns it as a string
%________________________________________________________________________________________________________________________

a = fopen([filename '.m']);
codeNumbers = fread(a);
theCode = char(codeNumbers);

end
