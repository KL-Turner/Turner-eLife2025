% this function was initially written to apply Kalman filter on the video
% tracking done by deeplabcut. 
% Now the main purpose of this code is to load the tracking data and do
% some housekeeping before saving the data. 
%% this function implements kalman filter
function KalmanFilterCorrection_OFT(filepath)
%% load the data
cd(filepath)
csvFiles = ls('*.csv'); % list the files in the folder
for csvList = 1:1:size(csvFiles,1)
        fileName = csvFiles(csvList,:);

        [csvData,bodyparts,rawcsv] = xlsread(fileName);
        %% Assign Data
        rmIdx(1) = 1;
        NumBodyParts = (size(bodyparts,2)-1)/3;
        for rr = 1:1:NumBodyParts
            rmIdx(rr+1)  = (rr*3)+1; 
        end
        BodyParts = bodyparts(2,:);
        TrackData = csvData;
        for il = 1:1:length(rmIdx)
            BodyParts(rmIdx(il)-(il-1)) = [];
            TrackData(:,rmIdx(il)-(il-1)) = [];
        end
        %% Run kalman filter to clean the data
%{
        % filter the non mouse points
%         for ObjectPoint = 4:1:NumBodyParts
%             RawData = TrackData(:,[(ObjectPoint*2)-1 ,ObjectPoint*2]);
%                 %% Create model parameters
%                 param.motionModel           = 'ConstantAcceleration';
%                 initialLocation             = RawData(1,:);
%                 param.initialEstimateError  = 1*ones(1,3); % use relatively small values
%                 param.motionNoise           = [25, 10, 1];
%                 param.measurementNoise      = 500;
% 
%                 %% create dummy dataset to overfit the model
%                 DummyData = RawData;
%                 rotationIdx = 200;
%                 for kk = 1:rotationIdx
%                 DummyData(length(RawData)+kk,:) = RawData(1,:);
%                 end
%                 DummyData = circshift(DummyData,rotationIdx);
% 
%                 %% train the model and correct the measurements
%                 for ll = 1: length(DummyData(:,1))  
%                     kalmanFilter = configureKalmanFilter(param.motionModel, ...
%                               initialLocation, param.initialEstimateError, ...
%                               param.motionNoise, param.measurementNoise);
%             %                 predict(kalmanFilter);
%                     trackedLocation(ll,:) = correct(kalmanFilter, DummyData(ll,:));
%                 end
% %                 FinalLocation = RawData;
%                 FinalLocation = trackedLocation(rotationIdx+1:length(DummyData(:,1)),:);
%                 CorrectTrackData(:,[(ObjectPoint*2)-1 ,ObjectPoint*2]) = FinalLocation;
%                 clear RawData DummyData FinalLocation
%         end
        %% no filter the mouse points
%         for ObjectPoint = 1:1:NumBodyParts %3
%                 CorrectTrackData(:,[(ObjectPoint*2)-1 ,ObjectPoint*2]) = TrackData(:,[(ObjectPoint*2)-1 ,ObjectPoint*2]);
%         end
%}
        %% Notes
        OFTData_cell{csvList} = TrackData;
        CSVData_cell{csvList} = rawcsv;
end
           
        TrackdataMat = cell2mat( OFTData_cell');
       
        extInd = strfind(fileName(1,:),'-'); 
     
        IdInd = strfind(fileName(1,:),'_'); 

        matfileFolder = ['..\..\..\Recalculated_Raw_Population\' fileName(1:IdInd (1)-1)]; 
        if ~isfolder(matfileFolder)
            mkdir(matfileFolder)
        end
        matfileName = [matfileFolder '\' fileName(1:extInd(1)-1) '.mat']; 
 
        OFTData.MouseID = fileName(1:IdInd(1)-1);
        OFTData.TrackData_cell = OFTData_cell;
        OFTData.TrackData = TrackdataMat;
        OFTData.bodyParts = BodyParts;
        OFTData.rawCSV = CSVData_cell;
        save(matfileName,'OFTData')

        % don't need to run this code now