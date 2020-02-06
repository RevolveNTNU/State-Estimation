function [] = SaveAsMat(analyze_data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dir_mat = uigetdir('', 'Select location to save .mat files');

channels = fieldnames(analyze_data);

for i = 1 : numel(channels)
    arr = analyze_data.(channels{i});
    save([dir_mat '\' channels{i}], 'arr');
end

end