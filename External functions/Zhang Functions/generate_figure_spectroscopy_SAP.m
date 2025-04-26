function generate_figure_spectroscopy_SAP()
close all; clear all; clc;

%% Data analysis
rootfolder = pwd; % the location of all small dataset

%% get the demographic info of all animals   
experiment.Blank = {
'T192';
'T205';
'T206';
'T208';
'T209';
'T211';
'T225';
};
experiment.SSP = {
'T200';
'T212';
'T213';
'T215';
'T216';
'T217';
'T218';
'T219';
};

%% set groups based on toxin injection
group.Blank.HR = [];
group.Blank.HbT = [];
group.Blank.HbD = [];
group.Blank.RestVar_HR = [];
group.Blank.RestVar_HbT = [];
group.Blank.RestVar_HbD = [];
group.SSP = group.Blank;

fields = fieldnames(experiment);

for a0 = 1:numel(fields)
    expCondition = fields{a0};
    for a1 = 1:numel(experiment.(expCondition))
        animalID = experiment.(expCondition){a1};
        disp(animalID);
        searchfolder = fullfile(rootfolder, 'Results', animalID);
        ind_targetFile = getfilenames(searchfolder, [expCondition,'-SAP', '*.mat']);
        if ~isempty(ind_targetFile)
            out.(animalID) = genFigure_individual_SAP(ind_targetFile);
        end
        group.(expCondition).HR = [group.(expCondition).HR;nanmean(out.(animalID).HR.LTA,1)];
        group.(expCondition).HbT = [group.(expCondition).HbT;nanmean(out.(animalID).HbT.PC.LTA,1)];
        group.(expCondition).HbD = [group.(expCondition).HbD;nanmean(out.(animalID).HbD.PC.LTA,1)];
        
        group.(expCondition).RestVar_HR = [group.(expCondition).RestVar_HR; nanmean(out.(animalID).RestVar.PC.HR)];
        group.(expCondition).RestVar_HbT = [group.(expCondition).RestVar_HbT; nanmean(out.(animalID).RestVar.PC.HbT)];
        group.(expCondition).RestVar_HbD = [group.(expCondition).RestVar_HbD; nanmean(out.(animalID).RestVar.PC.HbD)];
    end
end


% HR
figure_LTA(nanmean(group.Blank.HR,1),nanmean(group.SSP.HR,1),...
           nanstd(group.Blank.HR,[],1)/sqrt(size(group.Blank.HR,1)),...
           nanstd(group.SSP.HR,[],1)/sqrt(size(group.SSP.HR,1)));
legend({'Blank','SSP'});
ylabel('Heart rate (Hz)');

% oxygen       
figure_LTA(nanmean(group.Blank.HbD,1),nanmean(group.SSP.HbD,1),...
           nanstd(group.Blank.HbD,[],1)/sqrt(size(group.Blank.HbD,1)),...
           nanstd(group.SSP.HbD,[],1)/sqrt(size(group.SSP.HbD,1)));
legend({'Blank','SSP'});
ylabel('Oxygen, PC (uM)')


% HbT      
figure_LTA(nanmean(group.Blank.HbT,1),nanmean(group.SSP.HbT,1),...
           nanstd(group.Blank.HbT,[],1)/sqrt(size(group.Blank.HbT,1)),...
           nanstd(group.SSP.HbT,[],1)/sqrt(size(group.SSP.HbT,1)));
legend({'Blank','SSP'});
ylabel('HbT, PC (uM)')


%% HR/HbT/HbD variance, PC
x = [1.3 1.7];
figure('position',[100 100 800 400]);
h1 = subplot(131);
plot(x(1),group.Blank.RestVar_HR,...
    'marker','o','markersize',10,...
    'markeredgecolor','k','markerfacecolor',[0.5 0.5 0.5],'linestyle','none',...
    'color',[0.5 0.5 0.5]); % data points in FL/HL
hold on;
plot(x(2),group.SSP.RestVar_HR,...
    'marker','o','markersize',10,...
    'markeredgecolor','k','markerfacecolor','r','linestyle','none',...
    'color',[0.5 0.5 0.5]); % data points in FL/HL

bar(x(1), nanmean(group.Blank.RestVar_HR), 'BarWidth', 0.3,'EdgeColor','k', 'LineWidth', 1,'FaceColor','none');
bar(x(2), nanmean(group.SSP.RestVar_HR), 'BarWidth', 0.3,'EdgeColor','r', 'LineWidth', 1,'FaceColor','none');
% Using the following function to make sure the horizontal error bar the
% same size
errbar_QZ(x(1), nanmean(group.Blank.RestVar_HR), nanstd(group.Blank.RestVar_HR), 'k');
errbar_QZ(x(2), nanmean(group.SSP.RestVar_HR), nanstd(group.SSP.RestVar_HR), 'r');
axis([1 2 0 1]);
set(gca,'box','off','xtick', x, 'xticklabel',{'Blank','SSP'});
ylabel('HR SD, PC (Hz)');

h2 = subplot(132);
plot(x(1),group.Blank.RestVar_HbT,...
    'marker','o','markersize',10,...
    'markeredgecolor','k','markerfacecolor',[0.5 0.5 0.5],'linestyle','none',...
    'color',[0.5 0.5 0.5]); % data points in FL/HL
hold on;
plot(x(2),group.SSP.RestVar_HbT,...
    'marker','o','markersize',10,...
    'markeredgecolor','k','markerfacecolor','r','linestyle','none',...
    'color',[0.5 0.5 0.5]); % data points in FL/HL

bar(x(1), nanmean(group.Blank.RestVar_HbT), 'BarWidth', 0.3,'EdgeColor','k', 'LineWidth', 1,'FaceColor','none');
bar(x(2), nanmean(group.SSP.RestVar_HbT), 'BarWidth', 0.3,'EdgeColor','r', 'LineWidth', 1,'FaceColor','none');
% Using the following function to make sure the horizontal error bar the
% same size
errbar_QZ(x(1), nanmean(group.Blank.RestVar_HbT), nanstd(group.Blank.RestVar_HbT), 'k');
errbar_QZ(x(2), nanmean(group.SSP.RestVar_HbT), nanstd(group.SSP.RestVar_HbT), 'r');
axis([1 2 0 10]);
set(gca,'box','off','xtick', x, 'xticklabel',{'Blank','SSP'});
ylabel('HbT SD, PC (uM)');
 
% nanmean(group.Blank.RestVar_HbT)
% nanstd(group.Blank.RestVar_HbT) % 4.50 +/- 1.47
% 
% nanmean(group.SSP.RestVar_HbT)
% nanstd(group.SSP.RestVar_HbT) % 2.97 +/- 1.03
% 
% normal_distribution = adtest([group.Blank.RestVar_HbT;group.SSP.RestVar_HbT;]); % 0, normal
% equal_variance = vartest2(group.Blank.RestVar_HbT, group.SSP.RestVar_HbT); % 0, equal variance
% [h,p,ci,stats] = ttest2(group.Blank.RestVar_HbT(:),group.SSP.RestVar_HbT(:)); % P = 0.0338, t(13) = 2.3714


h3 = subplot(133);
plot(x(1),group.Blank.RestVar_HbD,...
    'marker','o','markersize',10,...
    'markeredgecolor','k','markerfacecolor',[0.5 0.5 0.5],'linestyle','none',...
    'color',[0.5 0.5 0.5]); % data points in FL/HL
hold on;
plot(x(2),group.SSP.RestVar_HbD,...
    'marker','o','markersize',10,...
    'markeredgecolor','k','markerfacecolor','r','linestyle','none',...
    'color',[0.5 0.5 0.5]); % data points in FL/HL

bar(x(1), nanmean(group.Blank.RestVar_HbD), 'BarWidth', 0.3,'EdgeColor','k', 'LineWidth', 1,'FaceColor','none');
bar(x(2), nanmean(group.SSP.RestVar_HbD), 'BarWidth', 0.3,'EdgeColor','r', 'LineWidth', 1,'FaceColor','none');
% Using the following function to make sure the horizontal error bar the
% same size
errbar_QZ(x(1), nanmean(group.Blank.RestVar_HbD), nanstd(group.Blank.RestVar_HbD), 'k');
errbar_QZ(x(2), nanmean(group.SSP.RestVar_HbD), nanstd(group.SSP.RestVar_HbD), 'r');
axis([1 2 0 20]);
set(gca,'box','off','xtick', x, 'xticklabel',{'Blank','SSP'});
ylabel('HbD SD, PC (uM)');

% nanmean(group.Blank.RestVar_HbD)
% nanstd(group.Blank.RestVar_HbD) % 11.98 +/- 2.25
% 
% nanmean(group.SSP.RestVar_HbD)
% nanstd(group.SSP.RestVar_HbD) % 10.21 +/- 1.15
% 
% normal_distribution = adtest([group.Blank.RestVar_HbD;group.SSP.RestVar_HbD;]); % 0, normal
% equal_variance = vartest2(group.Blank.RestVar_HbD, group.SSP.RestVar_HbD); % 0, equal variance
% [h,p,ci,stats] =
% ttest2(group.Blank.RestVar_HbD(:),group.SSP.RestVar_HbD(:)); % P = 0.0723, t(13) = 1.9558
