function new_data = Interpolate( data, timestamp_new )
%INTERPOLATE Summary of this function goes here
%   Detailed explanation goes here

new_data = interp1(data(:,1),data(:,2),timestamp_new);

end

