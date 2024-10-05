function [UniqueDays, dayIndex, dayID] = GetUniqueDays_2P(DateList)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpse: Takes a list of fileDates and determines how many unique individual days there are. 
%________________________________________________________________________________________________________________________
%
%   Inputs: A list of dates (YYMMDD)
%
%   Outputs: Number of unique days in strDay format as well as the index.
%
%   Last Revised: March 22nd, 2019
%________________________________________________________________________________________________________________________

if iscellstr(DateList)
    temp = cell2mat(DateList);
    DateList = temp;
end

filebreaks = strfind(DateList(1,:),'_');
if isempty(filebreaks)
    AllDates = DateList;
elseif or(length(filebreaks)==3,length(filebreaks)==4)
    AllDates = DateList(:,1:filebreaks(1)-1);
elseif length(filebreaks)==6
    date_ind = filebreaks(2)+1:filebreaks(3)-1;
    AllDates = DateList(:,date_ind);
else
    error('Format of the list of dates not recognized...')
end
All_days = mat2cell(AllDates,ones(1,size(AllDates,1)));
[UniqueDays, dayIndex, dayID] = unique(All_days);

end
