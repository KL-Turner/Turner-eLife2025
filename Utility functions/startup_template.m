function [] = startup_template()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% go to code directory - generate and add file paths of desired projects
cd C:\Users\kturner8\Documents
addpath(genpath('C:\Users\kturner8\Documents\Github\Data-Analysis'))
% addpath(genpath('C:\Users\kturner8\Documents\Github\Turner-eLife2020'))
% addpath(genpath('C:\Users\kturner8\Documents\Github\Turner-JNeurosci2023'))
% remove useless warnings
id = 'signal:filtfilt:ParseSOS';
warning('off',id)
clear id 
% add old style for figure toolbar because the new version sucks
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
% stop and enter debug upon errors
dbstop if error
% email settings
server = 'smtp.gmail.com';
mail = 'putemailhere@gmail.com';
password = 'putpasswordhere';
Apply prefs and props
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server',server);
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.port','587');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');