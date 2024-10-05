function [correctedTestData,FitParams] = CorrectKernelError_IOS_Manuscript2020(FitData,TestData,behavior)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner - originally written By Kyle Gheres
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: corrects the nonlinearity error from the HbT kernel
%________________________________________________________________________________________________________________________

if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
    fitPredArray = [];
    fitActualArray = [];
    % combine the fit data into a single array
    for aa = 1:length(FitData.mPred)
        fitPredArray = vertcat(fitPredArray,FitData.mPred{aa}'); %#ok<*AGROW>
        fitActualArray = vertcat(fitActualArray,FitData.mAct{aa}');
    end
    % combine the test data into a single array
    testPredArray = [];
    testActualArray = [];
    for bb = 1:length(TestData.mPred)
        testPredArray = vertcat(testPredArray,TestData.mPred{bb}');
        testActualArray = vertcat(testActualArray,TestData.mAct{bb}');
    end
elseif strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true
    % combine the fit data into a single array
    for aa = 1:size(FitData.mPred)
        fitPred(aa,:) = FitData.mPred{aa};
        fitActual(aa,:) = FitData.mAct{aa};
    end
    fitActualArray = reshape(fitActual,numel(fitActual),1);
    fitPredArray = reshape(fitPred,numel(fitPred),1);
    % combine the test data into a single array
    for bb = 1:size(TestData.mPred)
        testPred(bb,:) = TestData.mPred{bb};
        testActual(bb,:) = TestData.mAct{bb};
    end
    testActualArray = reshape(testActual,numel(testActual),1);
    testPredArray = reshape(testPred,numel(testPred),1);
end
% Correct non-linearity between predicted and actual HbT with respect to predicted HbT
[correctedTestPredArray,FitParams] = CorrectExpError_IOS_Manuscript2020(fitActualArray,fitPredArray,testActualArray,testPredArray);
% reshape data back to orignial form
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
    startPoint = 1;
    for a = 1:length(TestData.mPred)
        arrayLength = length(TestData.mPred{1,a});
        stopPoint = (startPoint + arrayLength) - 1;
        correctedTestData{a,1} = correctedTestPredArray(startPoint:stopPoint)';
        startPoint = startPoint + arrayLength;
    end 
elseif strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true
    correctedTestData = reshape(correctedTestPredArray,size(testPred,1),size(testPred,2));
end

end
