function [DataStruct,filtArray] = SelectConvolutionBehavioralEvents_HRF2020(DataStruct,behavior,hemisphere)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%
% Purpose: 
%________________________________________________________________________________________________________________________

BehaviorBufferTime = 4; % Seconds, 4 seconds allows CBV signals to stabalize
switch behavior
    case 'Rest'
        Criteria.Fieldname = {'PuffDistance','duration'};
        Criteria.Comparison = {'gt','gt'};
        Criteria.Value = {5,6 + BehaviorBufferTime};
        Criteria.Min_Duration = 5;
        Criteria.Min_PuffDist = 5;
    case 'Whisk'
        Criteria.Fieldname = {'duration','duration','restTime','PuffDistance'};
        Criteria.Comparison = {'lt','gt','gt','gt'};
        Criteria.Value = {3,0.5,BehaviorBufferTime,5};
        DataStruct = DataStruct.whisk;
    case 'Contra'
        Criteria.Fieldname = {'Name','whiskScore_Pre', 'movementScore_Pre','PuffDistance'};
        Criteria.Comparison = {'equal','lt','lt','gt'};
        Criteria.Value = {'Contra',0.1, 0.1, 5};
        % Criteria.Fieldname = {'Name'};
        % Criteria.Comparison = {'equal'};
        % Criteria.Value = {'Contra'};
        DataStruct = DataStruct.stim;
end

FName = Criteria.Fieldname;
Comp = Criteria.Comparison;
Val = Criteria.Value;

if length(FName)~=length(Comp)
    error(' ')
elseif length(FName)~=length(Val)
    error(' ')
end

if strcmp(behavior,'Contra')
    if strcmp(hemisphere,'LH') == true
        for a = 1:length(DataStruct.solenoidName)
            if strcmp(DataStruct.solenoidName{a,1},'LPadSol')
                DataStruct.Name{a,1} = 'Ipsi';
            elseif strcmp(DataStruct.solenoidName{a,1},'RPadSol')
                DataStruct.Name{a,1} = 'Contra';
            elseif strcmp(DataStruct.solenoidName{a,1},'AudSol')
                DataStruct.Name{a,1} = 'Control';
            end
        end
    elseif strcmp(hemisphere,'RH')
        for a = 1:length(DataStruct.solenoidName)
            if strcmp(DataStruct.solenoidName{a,1},'RPadSol')
                DataStruct.Name{a,1} = 'Ipsi';
            elseif strcmp(DataStruct.solenoidName{a,1},'LPadSol')
                DataStruct.Name{a,1} = 'Contra';
            elseif strcmp(DataStruct.solenoidName{a,1},'AudSol')
                DataStruct.Name{a,1} = 'Control';
            end
        end
    end
    DataStruct.duration = zeros(length(DataStruct.eventTime),1);
elseif strcmp(behavior,'Whisk')
    DataStruct.PuffDistance = DataStruct.puffDistance;
elseif strcmp(behavior,'Rest')
    DataStruct.PuffDistance = DataStruct.stimDistances;
    DataStruct.duration = DataStruct.durations;
    DataStruct.eventTime = DataStruct.eventTimes;
    DataStruct.samplingRate = DataStruct.samplingRate;
end

filtArray = true(size(DataStruct.data,1),1);
for FN = 1:length(FName)
    if ~isfield(DataStruct,FName{FN})
        error('Criteria field not found')
    end
    switch Comp{FN}
        case 'gt'
            if iscell(DataStruct.(FName{FN}))
                if ischar(DataStruct.(FName{FN}){1})
                    error(' ')
                else
                    IndFilt = false(size(filtArray));
                    for c = 1:length(DataStruct.(FName{FN}))
                        IndFilt(c) = all(gt(abs(DataStruct.(FName{FN}){c}), Val{FN}));
                    end
                end
            else
                IndFilt = gt(DataStruct.(FName{FN}), Val{FN});
            end
        case 'lt'
             if iscell(DataStruct.(FName{FN}))
                if ischar(DataStruct.(FName{FN}){1})
                    error(' ')
                else
                    IndFilt = false(size(filtArray));
                    for c = 1:length(DataStruct.(FName{FN}))
                        IndFilt(c) = all(lt(DataStruct.(FName{FN}){c}, Val{FN}));
                    end
                end
            else
                IndFilt = lt(DataStruct.(FName{FN}), Val{FN});
            end
        case 'equal'
            if iscell(DataStruct.(FName{FN}))
                IndFilt = strcmp(DataStruct.(FName{FN}),Val{FN});
            else
                IndFilt = DataStruct.(FName{FN}) == Val{FN};
            end
        otherwise
            error(' ')
    end
    % This is a patch for old data, only consider the pre-event movement
    if or(strcmp(FName{FN},'WhiskScore'),strcmp(FName{FN},'MoveScore'))
        IndFilt = IndFilt(:,1);
    end
    filtArray = and(filtArray,IndFilt);
end

end
