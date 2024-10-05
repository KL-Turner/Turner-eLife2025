function out = Cal_HbO_HbR(data, CMIdata, filter_data)
% Cal_HbO_HbR: calculate oxy- and deoxy-hemoglobin based on intensity
% changes from two different wavelengths
%
%Note:
% we use CMIdata just because in the paper from Hillman Lab, they only
% study resting state, therefore they can just do mean subtraction, and
% will get a stable baseline. But, for our data, we have running data, so
% we need to normalize the dataset before calculating hemoglobin
% concentration
%
% Inputs:
% data: raw image, speed data
% CMIdata: output from ComputeMeanIntensity.m, i.e., dR/R0 data
% filter_data: need to filter dataset? logical: true = filter below 1 Hz,
% false = raw dataset, 30 Hz
%
% Qingguang Zhang (qingguang.zhang@gmail.com)
% Last Modified: Oct. 5, 2019

%% Define pathlengths for each wavelength (cm)
% see following papers for more info.
% Ma et al., Wide-field optical mapping of neural activity and brain
% hemodynamics: considerations and novel approaches, Philos Trans R Soc
% Lond B Biol Sci, 371, 2016. [Ref #1]
%
% Ma et al., Resting-state hemodynamics are spatiotemporally coupled to 
% synchronized and symmetric neural activity in excitatory neurons, Proc Natl Acad 
% Sci U S A 113, 2016. [Ref #2]

% these parameters are from Ref #1
X_Blue = 0.0526691; % pathlength for (470nm) cm
X_Green = 0.0371713; % pathlength for  (530nm) cm, this is the pathlength of the green LEDs

% this parameter is from Ref #1
X_Red = 0.4535122; % pathlength for (660nm) cm

%% Define extinction coefficients (cm^-1*M^-1)
E_HbR_Blue = 16156.4; % 470 nm
E_HbO_Blue = 33209.2;

E_HbR_Red = 3226.56; % 660 nm
E_HbO_Red = 319.6;

% based on the data from https://omlc.org/spectra/hemoglobin/summary.html,
% the extinction coefficient for HbO and HbR should be different for 530 nm
% wavelength, but since we treat the 530nm as isobetic point, we set the
% extinction coefficient to be the same for HbO and HbR
% see Fig.2 in  Ma et al., Philos Trans R Soc Lond B Biol Sci, 371, 2016.

% E_HbO_Green = 39956.8; 
% E_HbR_Green = 39036.4; 
E_HbR_Green = 39500; 
E_HbO_Green = 39500;

%% Get intensity change of each pixel for each wavelength
if filter_data % use filtered data, i.e., < 1 Hz
    if data.Fr == 20 || data.Fr == 10
        CBV_Blue = [CMIdata.Blue.dRR0filtall{1};CMIdata.Blue.dRR0filtall{2};]+1;
        CBV_Red = [CMIdata.Red.dRR0filtall{1};CMIdata.Red.dRR0filtall{2};]+1;
        CBV_Green = [CMIdata.Green.dRR0filtall{1};CMIdata.Green.dRR0filtall{2}]+1;
    elseif data.Fr == 30
        CBV_Blue = [CMIdata.Blue.dRR0filtall{1};CMIdata.Blue.dRR0filtall{2};]+1;
        CBV_Red = [];
        CBV_Green = [CMIdata.Green.dRR0filtall{1};CMIdata.Green.dRR0filtall{2};]+1;
    end
else % for future correlation analysis between resp and HbD, use unfiltered data to get 30 Hz fluctuations
    if data.Fr == 20 || data.Fr == 10
        CBV_Blue = [CMIdata.Blue.dRR0all{1};CMIdata.Blue.dRR0all{2};]+1;
        CBV_Red = [CMIdata.Red.dRR0all{1};CMIdata.Red.dRR0all{2};]+1;
        CBV_Green = [CMIdata.Green.dRR0all{1};CMIdata.Green.dRR0all{2};]+1;
    elseif data.Fr == 30
        CBV_Blue = [CMIdata.Blue.dRR0all{1};CMIdata.Blue.dRR0all{2};]+1;
        CBV_Red = [];
        CBV_Green = [CMIdata.Green.dRR0all{1};CMIdata.Green.dRR0all{2};]+1;
    end
end

%% Calculate absorption coefficient for each pixel = -1/pathlength*ln(deltaR/R+1)
Mu_Blue = -1/X_Blue*log(CBV_Blue); % cm^-1, natural logarithm
Mu_Green = -1/X_Green*log(CBV_Green); % cm^-1, natural logarithm
Mu_Red = -1/X_Red*log(CBV_Red); % cm^-1, natural logarithm

%% Get intensity change for each ROI for each wavelength
% preallocate for speed
CBV_Blue_ROI = cell(12,1);
CBV_Red_ROI = cell(12,1);
CBV_Green_ROI = cell(12,1);

if filter_data % use filtered data
    if data.Fr == 20 || data.Fr == 10
        for i = 1:numel(CMIdata.Blue.dRR0meanall)
            CBV_Blue_ROI{i} = CMIdata.Blue.dRR0meanall{i}+1;
            CBV_Red_ROI{i} = CMIdata.Red.dRR0meanall{i}+1;
            CBV_Green_ROI{i} = CMIdata.Green.dRR0meanall{i}+1;
        end
    elseif data.Fr == 30
        for i = 1:numel(CMIdata.Blue.dRR0meanall)
            CBV_Blue_ROI{i} = CMIdata.Blue.dRR0meanall{i}+1;
            CBV_Red_ROI{i} = [];
            CBV_Green_ROI{i} = CMIdata.Green.dRR0meanall{i}+1;
        end
    end
else % use unfiltered data
    if data.Fr == 20 || data.Fr == 10
        for i = 1:numel(CMIdata.Blue.dRR0meanall)
            CBV_Blue_ROI{i} = mean(CMIdata.Blue.dRR0all{i})+1;
            CBV_Red_ROI{i} = mean(CMIdata.Red.dRR0all{i})+1;
            CBV_Green_ROI{i} = mean(CMIdata.Green.dRR0all{i})+1;
        end
    elseif data.Fr == 30
        for i = 1:numel(CMIdata.Blue.dRR0meanall)
            CBV_Blue_ROI{i} = mean(CMIdata.Blue.dRR0all{i})+1;
            CBV_Red_ROI{i} = [];
            CBV_Green_ROI{i} = mean(CMIdata.Green.dRR0all{i})+1;
        end
    end
end

%% Calculate absorption coefficient for each ROI = -1/pathlength*ln(deltaR/R+1)
% preallocate for speed
Mu_Blue_ROI = cell(12,1);
Mu_Red_ROI = cell(12,1);
Mu_Green_ROI = cell(12,1);

for i = 1:numel(CBV_Blue_ROI)
    Mu_Blue_ROI{i} = -1/X_Blue*log(CBV_Blue_ROI{i}); % cm^-1, natural logarithm
    Mu_Green_ROI{i} = -1/X_Green*log(CBV_Green_ROI{i}); % cm^-1, natural logarithm
    Mu_Red_ROI{i} = -1/X_Red*log(CBV_Red_ROI{i}); % cm^-1, natural logarithm
end

%% Calculate concentration of HbR & HbO for each pixel
% Calculate concentrations (uM) using blue and green light
C_HbR.BG = (E_HbO_Green*Mu_Blue - E_HbO_Blue*Mu_Green)/...
    (E_HbO_Green*E_HbR_Blue - E_HbO_Blue*E_HbR_Green)*1e6; % in uM
C_HbO.BG = (E_HbR_Green*Mu_Blue - E_HbR_Blue*Mu_Green)/...
    (E_HbR_Green*E_HbO_Blue - E_HbR_Blue*E_HbO_Green)*1e6; % in uM
D_HbO_HbR.BG = C_HbO.BG - C_HbR.BG; % difference between HbO and HbR

if ~isempty(CBV_Red)
    % Calculate concentrations (uM) using green and red light
    C_HbR.RG = (E_HbO_Green*Mu_Red - E_HbO_Red*Mu_Green)/...
        (E_HbO_Green*E_HbR_Red - E_HbO_Red*E_HbR_Green)*1e6; % in uM
    C_HbO.RG = (E_HbR_Green*Mu_Red - E_HbR_Red*Mu_Green)/...
        (E_HbR_Green*E_HbO_Red - E_HbR_Red*E_HbO_Green)*1e6; % in uM
    D_HbO_HbR.RG = C_HbO.RG - C_HbR.RG; % difference between HbO and HbR

    % Calculate concentrations (uM) using blue and red light
    C_HbR.BR = (E_HbO_Blue*Mu_Red - E_HbO_Red*Mu_Blue)/...
        (E_HbO_Blue*E_HbR_Red - E_HbO_Red*E_HbR_Blue)*1e6; % in uM
    C_HbO.BR = (E_HbR_Blue*Mu_Red - E_HbR_Red*Mu_Blue)/...
        (E_HbR_Blue*E_HbO_Red - E_HbR_Red*E_HbO_Blue)*1e6; % in uM  
    D_HbO_HbR.BR = C_HbO.BR - C_HbR.BR; % difference between HbO and HbR
else
    % Calculate concentrations (uM) using green and red light
    C_HbR.RG = [];
    C_HbO.RG = [];
    D_HbO_HbR.RG = []; % difference between HbO and HbR

    % Calculate concentrations (uM) using blue and red light
    C_HbR.BR = [];
    C_HbO.BR = [];
    D_HbO_HbR.BR = []; % difference between HbO and HbR
end

% HbO and HbR for each pixel
out.C_HbR = C_HbR;
out.C_HbO = C_HbO;
out.D_HbO_HbR = D_HbO_HbR;

%% Get concentration of HbR & HbO for each ROI
% preallocate for speed
C_HbR_ROI = cell(12,1);
C_HbO_ROI = cell(12,1);
D_HbO_HbR_ROI = cell(12,1);

for i = 1:numel(Mu_Blue_ROI)
    % Calculate concentrations (uM) using blue and green light
    C_HbR_ROI{i}.BG = (E_HbO_Green*Mu_Blue_ROI{i} - E_HbO_Blue*Mu_Green_ROI{i})/...
        (E_HbO_Green*E_HbR_Blue - E_HbO_Blue*E_HbR_Green)*1e6; % in uM
    C_HbO_ROI{i}.BG = (E_HbR_Green*Mu_Blue_ROI{i} - E_HbR_Blue*Mu_Green_ROI{i})/...
        (E_HbR_Green*E_HbO_Blue - E_HbR_Blue*E_HbO_Green)*1e6; % in uM
    D_HbO_HbR_ROI{i}.BG = C_HbO_ROI{i}.BG - C_HbR_ROI{i}.BG; % difference between HbO and HbR
    
    if ~isempty(CBV_Red)
        % Calculate concentrations (uM) using green and red light
        C_HbR_ROI{i}.RG = (E_HbO_Green*Mu_Red_ROI{i} - E_HbO_Red*Mu_Green_ROI{i})/...
            (E_HbO_Green*E_HbR_Red - E_HbO_Red*E_HbR_Green)*1e6; % in uM
        C_HbO_ROI{i}.RG = (E_HbR_Green*Mu_Red_ROI{i} - E_HbR_Red*Mu_Green_ROI{i})/...
            (E_HbR_Green*E_HbO_Red - E_HbR_Red*E_HbO_Green)*1e6; % in uM
        D_HbO_HbR_ROI{i}.RG = C_HbO_ROI{i}.RG - C_HbR_ROI{i}.RG; % difference between HbO and HbR
       
        % Calculate concentrations (uM) using blue and red light
        C_HbR_ROI{i}.BR = (E_HbO_Blue*Mu_Red_ROI{i} - E_HbO_Red*Mu_Blue_ROI{i})/...
            (E_HbO_Blue*E_HbR_Red - E_HbO_Red*E_HbR_Blue)*1e6; % in uM
        C_HbO_ROI{i}.BR = (E_HbR_Blue*Mu_Red_ROI{i} - E_HbR_Red*Mu_Blue_ROI{i})/...
            (E_HbR_Blue*E_HbO_Red - E_HbR_Red*E_HbO_Blue)*1e6; % in uM
        D_HbO_HbR_ROI{i}.BR = C_HbO_ROI{i}.BR - C_HbR_ROI{i}.BR; % difference between HbO and HbR
    else
        % Calculate concentrations (uM) using green and red light
        C_HbR_ROI{i}.RG = [];
        C_HbO_ROI{i}.RG = [];
        D_HbO_HbR_ROI{i}.RG = [];
    
        % Calculate concentrations (uM) using blue and red light
        C_HbR_ROI{i}.BR = [];
        C_HbO_ROI{i}.BR = [];
        D_HbO_HbR_ROI{i}.BR = [];
    end
end

out.C_HbR_ROI = C_HbR_ROI;
out.C_HbO_ROI = C_HbO_ROI;
out.D_HbO_HbR_ROI = D_HbO_HbR_ROI;

end