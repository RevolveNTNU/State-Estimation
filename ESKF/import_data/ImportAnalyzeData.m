function [analyzeData] = ImportAnalyzeData(folder)

if nargin < 1
    folder =  uigetdir();
end

% Find all csv files
dataFiles = dir(strcat(folder, '\*.csv'));
fileNames = {dataFiles(:).name};

% The channel used to interpolate all other channels and get the time from

% 2019
masterChannelFile = 'vcu_INS_ax.csv';

% % 2018
% masterChannelFile = 'SBS_F1_Steering_Angle.csv';

% masterChannelFile = 'FL_encoder_angle.csv';



%Check if the masterchannel exists
idx = find(strcmp(fileNames, masterChannelFile));
if not(isnumeric(idx))
    error('Cant find masterChannel');
end

% Remove tha masterchannel from the list of regular channels
fileNames(:,idx) = [];

analyzeData = struct;


masterChannelData = importChannel(strcat(folder, '\', masterChannelFile));
analyzeData.('time') = masterChannelData(:,1);
analyzeData.(masterChannelFile(1:end-4)) = masterChannelData(:,2); % (1:end-4) trims avay file extension to get varName


for i =1:length(fileNames)
    fileName = fileNames{i};
    fileName
    fileData = importChannel(strcat(folder,'\', fileName));
    interpolatedValues =  Interpolate(fileData, analyzeData.('time'));
    varName = fileName(1:end-4); % (1:end-4) trims avay file extension to get varName
    varName = strrep(varName, '.', '_');
    analyzeData.(varName) = interpolatedValues;
end