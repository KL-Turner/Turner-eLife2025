function [] = CreateTrainingDataSet_2P_eLife2025(mergedDataFileIDs,RestingBaselines,baselineType)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
for a = 1:size(mergedDataFileIDs,1)
    mergedDataFileID = mergedDataFileIDs(a,:);
    trainingDataFileID = [mergedDataFileID(1:end - 14) 'TrainingData.mat'];
    if ~exist(trainingDataFileID,'file')
        disp(['Loading ' mergedDataFileID ' for manual sleep scoring.' ]); disp(' ')
        load(mergedDataFileID)
        saveFigs = 'n';
        [figHandle] = GenerateSingleFigures_2P_eLife2025(mergedDataFileID,baselineType,saveFigs,RestingBaselines);
        trialDuration = MergedData.notes.trialDuration_Sec;
        numBins = trialDuration/5;
        behavioralState = cell(180,1);
        for b = 1:numBins
            global buttonState %#ok<TLEV>
            buttonState = 0;
            xStartVal = (b*5) - 4;
            xEndVal = b*5;
            xInds = xStartVal:1:xEndVal;
            figHandle = gcf;
            subplot(6,1,4)
            hold on
            leftEdge4 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
            rightEdge4 = xline(xInds(5),'color',colors('electric purple'),'LineWidth',2);
            if b <= 60
                xlim([1,300])
            elseif b >= 61 && b <= 120
                xlim([300,600])
            elseif b >= 121 && b <= 180
                xlim([600,900])
            end
            [updatedGUI] = SelectBehavioralStateGUI;
            while buttonState == 0
                drawnow()
                if buttonState == 1
                    guiResults = guidata(updatedGUI);
                    if guiResults.togglebutton1.Value == true
                        behavioralState{b,1} = 'Not Sleep';
                    elseif guiResults.togglebutton2.Value == true
                        behavioralState{b,1} = 'NREM Sleep';
                    elseif guiResults.togglebutton3.Value == true
                        behavioralState{b,1} = 'REM Sleep';
                    else
                        disp('No button pressed'); disp(' ')
                        keyboard
                    end
                    close(updatedGUI)
                    break;
                end
                ...
            end
        delete(leftEdge4)
        delete(rightEdge4)
        end
        close(figHandle)
        paramsTable.behavState = behavioralState;
        trainingTable = paramsTable;
        save(trainingDataFileID,'trainingTable')
    else
        disp([trainingDataFileID ' already exists. Continuing...']); disp(' ')
    end
end

end
