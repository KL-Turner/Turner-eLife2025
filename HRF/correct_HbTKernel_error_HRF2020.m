function [FitParams,CorrectedHbT,TestActual,TestPred]=correct_HbTKernel_error(FitData,TestData)

%% Get Half the data to fit
splitCount=1;
for trialnum=1:2:length(SampleData.actualHbT)
    FitActual(splitCount,:)=SampleData.actualHbT{trialnum};
    FitPred(splitCount,:)=SampleData.predictedHbT{trialnum};
    splitCount=splitCount+1;
end
FitActual=reshape(FitActual,numel(FitActual),1);
FitPred=reshape(FitPred,numel(FitPred),1);
%% Get other half of data to test correction
splitCount=1;
for trialnum=2:2:length(SampleData.actualHbT)
    TestActual(splitCount,:)=SampleData.actualHbT{trialnum};
    TestPred(splitCount,:)=SampleData.predictedHbT{trialnum};
    splitCount=splitCount+1;
end
TestActual=reshape(TestActual,numel(TestActual),1);
TestPred=reshape(TestPred,numel(TestPred),1);

for a = 1:size(SampleData.actualHbT)
    FitActual(a,:) = SampleData.actualHbT{a};
    FitPred(a,:) = SampleData.predictedHbT{a};
    TestActual(a,:) = SampleData.actualHbT{a};
    TestPred(a,:) = SampleData.predictedHbT{a};
end
FitActual=reshape(FitActual,numel(FitActual),1);
FitPred=reshape(FitPred,numel(FitPred),1);
TestActual=reshape(TestActual,numel(TestActual),1);
TestPred=reshape(TestPred,numel(TestPred),1);

%% Correct non-linearity between predicted and actual HbT with respect to predicted HbT
[FitParams,CorrectedHbT]=correctExpError(FitActual,FitPred,TestActual,TestPred);
CorrectedHbT = reshape(CorrectedHbT,size(SampleData.actualHbT,1),size(SampleData.actualHbT{1,1},2));
close all

% %% Correct Linear error in prediction kernel data with respect to known IOS Data
% [FitParams,CorrectedHbT]=correctLinearError(FitActual,FitPred,TestActual,TestPred);
% 
% %% Correct error in prediction kernel data with respect to predicted HbT
% [FitParams,CorrectedHbT]=correctExpError_Predicted(FitActual,FitPred,TestActual,TestPred);
