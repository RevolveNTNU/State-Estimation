function [analyzeData] = ImportRaw(folder)
if nargin < 1
    folder =  uigetdir();
end

% Find all csv files
dataFiles = dir(strcat(folder, '\*.csv'));
fileNames = {dataFiles(:).name};

analyzeData = struct;



for i =1:length(fileNames)
    fileName = fileNames{i};
%     fileName
    fileData = importChannel(strcat(folder,'\', fileName));
    fileData(:,1) = fileData(:,1) * 1e-6;
    varName = fileName(1:end-4); % (1:end-4) trims avay file extension to get varName
    varName = strrep(varName, '.', '_');
    analyzeData.(varName) = fileData;
end
end

