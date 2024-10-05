function [CorrectedHbT,FitParams] = CorrectExpError_IOS_Manuscript2020(FitActual,FitPred,TestActual,TestPred)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner - originally written By Kyle Gheres
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
% 
% Purpose: This function takes HbT data taken from IOS experiments and paired predicted HbT data derived from a convolution
%          of an impulse response kernel and neural activity collected simultaneously with known IOS data, fits a two 
%          exponent function (A*exp(B*x(t))+C*exp(D*x(t)) to describe the relationship, and removes this trend from the
%          predicted IOS data. 
%________________________________________________________________________________________________________________________

close all;

%% Visualize Real vs Predicted data by trial
figure(89); hold on
scatter(FitPred,FitActual);
xlabel('Predicted \DeltaHbT (\muM)');
ylabel('Actual \DeltaHbT (\muM)');
title('Actual vs Predicted \Deltahemoglobin');

%% Visualize Residual vs Predicted data
figure(88);
scatter(FitActual,(FitPred-FitActual));
xlabel('Actual \DeltaHbT (\muM)');
ylabel('Residual (Predicted \DeltaHbT- Actual \DeltaHbT) (\muM)');
title('Residual vs Predicted \Deltahemoglobin');

%% Plot 2D histogram of Real vs Predicted
BinResolution=0.5; %bin size for histogram uM.
XEdges=-200:BinResolution:200;
YEdges=-200:BinResolution:200;

figure(101);
FitHist=histogram2(FitPred,FitActual,XEdges,YEdges);
xlabel('Predicted \DeltaHbT (\muM)');
ylabel('Actual \DeltaHbT (\muM)');
zlabel('Bin Counts');

%% Normalize bin counts for each Real HbT column by total number of counts in each column
theCounts=FitHist.BinCounts;
ColCounts=sum(theCounts,2);
TotCounts=sum(ColCounts);
NormMat=repmat(ColCounts,1,size(theCounts,2));
NormCounts=(theCounts./NormMat)*100;
NormCounts(isnan(NormCounts))=0;
TotBins=(theCounts/TotCounts)*100;
TotBins(isnan(TotBins))=0;

figure(102);
NormHist=imagesc(XEdges,YEdges,NormCounts');
caxis([0 5]);
colorbar('eastoutside')
axis xy;
xlabel('Predicted HbT (\muM)');
ylabel('Actual HbT (\muM)');
title('Column normalized histogram of Actual vs Predicted hemoglobin');

figure(99);
NormHist=imagesc(XEdges,YEdges,TotBins');
caxis([0 0.3]);
colorbar('eastoutside')
axis xy;
xlabel('Predicted HbT (\muM)');
ylabel('Actual HbT (\muM)');
title('Normalized histogram of Actual vs Predicted hemoglobin');

%% Find the bin with the highest number of counts for each real HbT bin
for colNum=1:size(NormCounts,1)
    maxVal=max(NormCounts(colNum,:));
    if maxVal~=0
        MaxBinLocs=find(NormCounts(colNum,:)==maxVal);
        PeakLocs(colNum)=round(median(MaxBinLocs),0);
    else
        PeakLocs(colNum)=NaN;
    end
end
PredictedPeaks(1:length(XEdges))=0;
PredictedPeaks(~isnan(PeakLocs))=YEdges(PeakLocs(~isnan(PeakLocs)));

figure(103); hold on;
imagesc(XEdges,YEdges,NormCounts');
caxis([0 5]);
colorbar('eastoutside')
axis xy;
xlabel('Predicted \DeltaHbT (\muM)');
ylabel('Observed \DeltaHbT (\muM)');
title('Normalized histogram of Actual vs Predicted hemoglobin');
plot(XEdges,PredictedPeaks,'r','LineWidth',2);
xlim([-200 200]);
ylim([-200 200]);
legend('Peak count at predicted \DeltaHbT (\muM)');


%% Visualize histogram of Residuals vs predicted data
figure(172);
ResHist=histogram2(FitActual,(FitPred-FitActual),XEdges,YEdges);
xlabel('Actual \DeltaHbT (\muM)');
ylabel('Residual \DeltaHbT (\muM)');
zlabel('Bin Counts');

ResCounts=ResHist.BinCounts;
ColCounts=sum(ResCounts,2);
NormResMat=repmat(ColCounts,1,size(theCounts,2));
NormResCounts=(ResCounts./NormResMat)*100; 
NormResCounts(isnan(NormResCounts))=0;

figure(173);
NormResHist=imagesc(XEdges,YEdges,NormResCounts');
caxis([0 5]);
colorbar('eastoutside')
axis xy;
xlabel('Actual \DeltaHbT (\muM)');
ylabel('Residual \DeltaHbT (\muM)');
title('Column normalized histogram of Residual vs Actual \DeltaHbT');

%% Fit Predicted Peaks with Exponential
EmptyBins=find(isnan(PeakLocs)==1);
StartPoint=EmptyBins(find(diff(EmptyBins)>5,1,'first'))+1;
EndPoint=EmptyBins(find(diff(EmptyBins)>5,1,'first')+1)-1;
theEqn='a*exp(b*x)+c*exp(d*x)+(e*x)';
StartPts=[20 0.005 -25 -0.05 1];

[curve,goodness,output]=fit(XEdges(StartPoint:EndPoint)',PredictedPeaks(StartPoint:EndPoint)','exp2','Robust','LAR');
FitParams.ActualvPredicted.fitFormula=formula(curve);
FitParams.ActualvPredicted.fitCoeffs=coeffvalues(curve);
FitParams.ActualvPredicted.fitNames=coeffnames(curve);
FitParams.ActualvPredicted.FitStats=goodness;
FitParams.ActualvPredicted.FitInfo=output;

[curve,goodness,output]=fit(FitPred,FitActual,'exp2','Robust','LAR');
FitParams.AllPoints.fitFormula=formula(curve);
FitParams.AllPoints.fitCoeffs=coeffvalues(curve);
FitParams.AllPoints.fitNames=coeffnames(curve);
FitParams.AllPoints.FitStats=goodness;
FitParams.AllPoints.FitInfo=output;

[curve,goodness,output]=fit(XEdges(StartPoint:EndPoint)',PredictedPeaks(StartPoint:EndPoint)',theEqn,'Start',StartPts,'Robust','LAR');
FitParams.CustomFitPeaks.fitFormula=formula(curve);
FitParams.CustomFitPeaks.fitCoeffs=coeffvalues(curve);
FitParams.CustomFitPeaks.fitNames=coeffnames(curve);
FitParams.CustomFitPeaks.FitStats=goodness;
FitParams.CustomFitPeaks.FitInfo=output;

[curve,goodness,output]=fit(FitPred,FitActual,theEqn,'Start',StartPts,'Robust','LAR');
FitParams.CustomFitPoints.fitFormula=formula(curve);
FitParams.CustomFitPoints.fitCoeffs=coeffvalues(curve);
FitParams.CustomFitPoints.fitNames=coeffnames(curve);
FitParams.CustomFitPoints.FitStats=goodness;
FitParams.CustomFitPoints.FitInfo=output;

%% Correct Kernel Predicted Data
CorrectedHbT=FitParams.ActualvPredicted.fitCoeffs(1)*exp(FitParams.ActualvPredicted.fitCoeffs(2)*TestPred)+FitParams.ActualvPredicted.fitCoeffs(3)*exp(FitParams.ActualvPredicted.fitCoeffs(4)*TestPred);

HistPred=FitParams.ActualvPredicted.fitCoeffs(1)*exp(FitParams.ActualvPredicted.fitCoeffs(2)*XEdges)+FitParams.ActualvPredicted.fitCoeffs(3)*exp(FitParams.ActualvPredicted.fitCoeffs(4)*XEdges);

ScatterPred=FitParams.ActualvPredicted.fitCoeffs(1)*exp(FitParams.ActualvPredicted.fitCoeffs(2)*(-50:0.25:200))+FitParams.ActualvPredicted.fitCoeffs(3)*exp(FitParams.ActualvPredicted.fitCoeffs(4)*(-50:0.25:200));

CorrectedHbT2=FitParams.AllPoints.fitCoeffs(1)*exp(FitParams.AllPoints.fitCoeffs(2)*TestPred)+FitParams.AllPoints.fitCoeffs(3)*exp(FitParams.AllPoints.fitCoeffs(4)*TestPred);

HistPred2=FitParams.AllPoints.fitCoeffs(1)*exp(FitParams.AllPoints.fitCoeffs(2)*XEdges)+FitParams.AllPoints.fitCoeffs(3)*exp(FitParams.AllPoints.fitCoeffs(4)*XEdges);

ScatterPred2=FitParams.AllPoints.fitCoeffs(1)*exp(FitParams.AllPoints.fitCoeffs(2)*(-50:0.25:200))+FitParams.AllPoints.fitCoeffs(3)*exp(FitParams.AllPoints.fitCoeffs(4)*(-50:0.25:200));

CorrectedHbT3=FitParams.CustomFitPoints.fitCoeffs(1)*exp(FitParams.CustomFitPoints.fitCoeffs(2)*TestPred)+FitParams.CustomFitPoints.fitCoeffs(3)*exp(FitParams.CustomFitPoints.fitCoeffs(4)*TestPred)...
    +(FitParams.CustomFitPoints.fitCoeffs(5)*TestPred);%+FitParams.CustomFitPoints.fitCoeffs(6);

HistPred3=FitParams.CustomFitPoints.fitCoeffs(1)*exp(FitParams.CustomFitPoints.fitCoeffs(2)*XEdges)+FitParams.CustomFitPoints.fitCoeffs(3)*exp(FitParams.CustomFitPoints.fitCoeffs(4)*XEdges)...
    +(FitParams.CustomFitPoints.fitCoeffs(5)*XEdges);%+FitParams.CustomFitPoints.fitCoeffs(6);

ScatterPred3=FitParams.CustomFitPoints.fitCoeffs(1)*exp(FitParams.CustomFitPoints.fitCoeffs(2)*(-50:0.25:200))+FitParams.CustomFitPoints.fitCoeffs(3)*exp(FitParams.CustomFitPoints.fitCoeffs(4)*(-50:0.25:200))...
    +(FitParams.CustomFitPoints.fitCoeffs(5)*(-50:0.25:200));%+FitParams.CustomFitPoints.fitCoeffs(6);

CorrectedHbT4=FitParams.CustomFitPeaks.fitCoeffs(1)*exp(FitParams.CustomFitPeaks.fitCoeffs(2)*TestPred)+FitParams.CustomFitPeaks.fitCoeffs(3)*exp(FitParams.CustomFitPeaks.fitCoeffs(4)*TestPred)...
    +(FitParams.CustomFitPeaks.fitCoeffs(5)*TestPred);%+FitParams.CustomFitPoints.fitCoeffs(6);

HistPred4=FitParams.CustomFitPeaks.fitCoeffs(1)*exp(FitParams.CustomFitPeaks.fitCoeffs(2)*XEdges)+FitParams.CustomFitPeaks.fitCoeffs(3)*exp(FitParams.CustomFitPeaks.fitCoeffs(4)*XEdges)...
    +(FitParams.CustomFitPeaks.fitCoeffs(5)*XEdges);%+FitParams.CustomFitPoints.fitCoeffs(6);

ScatterPred4=FitParams.CustomFitPeaks.fitCoeffs(1)*exp(FitParams.CustomFitPeaks.fitCoeffs(2)*(-50:0.25:200))+FitParams.CustomFitPeaks.fitCoeffs(3)*exp(FitParams.CustomFitPeaks.fitCoeffs(4)*(-50:0.25:200))...
    +(FitParams.CustomFitPeaks.fitCoeffs(5)*(-50:0.25:200));%+FitParams.CustomFitPoints.fitCoeffs(6);

figure(134); hold on;
imagesc(XEdges,YEdges,NormCounts');
caxis([0 5]);
colorbar('eastoutside')
axis xy;
xlabel('Predicted \DeltaHbT (\muM)');
ylabel('Observed \DeltaHbT (\muM)');
title('Normalized histogram of Actual vs Predicted hemoglobin');
plot(XEdges,HistPred,'r','LineWidth',2);
plot(XEdges,HistPred2,'m','LineWidth',2);
plot(XEdges,HistPred3,'c','LineWidth',2);
plot(XEdges,HistPred4,'g','LineWidth',2);
xlim([-50 200]);
ylim([-100 100]);
legend({'Corrected \DeltaHbT (\muM) peak counts','Corrected \DeltaHbT (\muM) all data',...
    'Corrected \DeltaHbT (\muM) w/linear component points','Corrected \DeltaHbT (\muM) w/linear component peaks'});

figure(139); hold on
scatter(FitPred,FitActual);
plot((-50:0.25:200),ScatterPred,'r','LineWidth',2);
plot((-50:0.25:200),ScatterPred2,'m','LineWidth',2);
plot((-50:0.25:200),ScatterPred3,'c','LineWidth',2);
plot((-50:0.25:200),ScatterPred4,'g','LineWidth',2);
xlabel('Predicted \DeltaHbT (\muM)');
ylabel('Actual \DeltaHbT (\muM)');
title('Actual \DeltaHbT vs Predicted \DeltaHbT used to generate fit');
xlim([-25 200]);

figure(140); hold on
scatter(TestPred,TestActual);
plot((-50:0.25:200),ScatterPred,'r','LineWidth',2);
plot((-50:0.25:200),ScatterPred2,'m','LineWidth',2);
plot((-50:0.25:200),ScatterPred3,'c','LineWidth',2);
plot((-50:0.25:200),ScatterPred4,'g','LineWidth',2);
xlabel('Predicted \DeltaHbT (\muM)');
ylabel('Actual \DeltaHbT (\muM)');
title('Actual \DeltaHbT vs Predicted \DeltaHbT to test fit');
xlim([-25 200]);

figure(50); hold on
scatter(CorrectedHbT,TestActual);
xlabel('Corrected \DeltaHbT (\muM)');
ylabel('Actual \DeltaHbT (\muM)');
title('Corrected \DeltaHbT vs Real \DeltaHbT ');

figure(51); hold on
scatter(CorrectedHbT2,TestActual);
xlabel('Corrected \DeltaHbT (\muM)');
ylabel('Actual \DeltaHbT (\muM)');
title('Corrected \DeltaHbT vs Real \DeltaHbT ');

figure(53); hold on
scatter(CorrectedHbT3,TestActual);
xlabel('Corrected \DeltaHbT (\muM)');
ylabel('Actual \DeltaHbT (\muM)');
title('Corrected \DeltaHbT vs Real \DeltaHbT ');

figure(54); hold on
scatter(CorrectedHbT4,TestActual);
xlabel('Corrected \DeltaHbT (\muM)');
ylabel('Actual \DeltaHbT (\muM)');
title('Corrected \DeltaHbT vs Real \DeltaHbT ');

figure(52); hold on
scatter(TestActual,(CorrectedHbT-TestActual));
xlabel('Actual \DeltaHbT (\muM)');
ylabel('Residual \DeltaHbT (\muM)');
title('Real vs Residuals of Corrected \DeltaHbT');

figure(55); hold on
scatter(TestActual,(CorrectedHbT3-TestActual));
xlabel('Actual \DeltaHbT (\muM)');
ylabel('Residual \DeltaHbT (\muM)');
title('Real vs Residuals of Corrected \DeltaHbT');

figure(56); hold on
scatter(TestActual,(CorrectedHbT4-TestActual));
xlabel('Actual \DeltaHbT (\muM)');
ylabel('Residual \DeltaHbT (\muM)');
title('Real vs Residuals of Corrected \DeltaHbT');

%% Normalized Histogram of Test Data
figure(144);
TestHist=histogram2(TestPred,TestActual,XEdges,YEdges);
xlabel('Predicted \DeltaHbT (\muM)');
ylabel('Actual \DeltaHbT (\muM)');
zlabel('Bin Counts');

%% Normalize bin counts for each Real HbT column by total number of counts in each column
theCounts=TestHist.BinCounts;
ColCounts=sum(theCounts,2);
NormMat=repmat(ColCounts,1,size(theCounts,2));
NormCounts=(theCounts./NormMat)*100;
NormCounts(isnan(NormCounts))=0;

figure(145); hold on;
imagesc(XEdges,YEdges,NormCounts');
caxis([0 5]);
colorbar('eastoutside')
axis xy;
xlabel('Predicted \DeltaHbT (\muM)');
ylabel('Observed \DeltaHbT (\muM)');
title('Normalized histogram of Actual vs Predicted hemoglobin');
plot(XEdges,HistPred,'r','LineWidth',2);
plot(XEdges,HistPred2,'m','LineWidth',2);
xlim([-50 200]);
ylim([-100 100]);

end
