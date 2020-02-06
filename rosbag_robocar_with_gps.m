clear all;
clc;
close all;

bag = rosbag('C:\Users\tonja\OneDrive\Documents\rosbags\robocar_dataset_20180502.bag');
bag.AvailableTopics;
bag.MessageList;
start = bag.StartTime;
imubagselect = select(bag,'Topic','/imu/data');

% ts = timeseries(bagselect, 'LinearAcceleration.X', 'LinearAcceleration.Y', 'LinearAcceleration.Z' );
% ts.Data;

msgs = readMessages(imubagselect,'DataFormat','struct');
% msgs{3};
ax = cellfun(@(m) double(m.LinearAcceleration.X), msgs);
ay = cellfun(@(m) double(m.LinearAcceleration.Y), msgs);
az = cellfun(@(m) double(m.LinearAcceleration.Z), msgs);
 
omegax = cellfun(@(m) double(m.AngularVelocity.X), msgs);
omegay = cellfun(@(m) double(m.AngularVelocity.Y), msgs);
omegaz = cellfun(@(m) double(m.AngularVelocity.Z), msgs);
 

gpsbagselect = select(bag,'Topic','/fix');
msgs = readMessages(gpsbagselect,'DataFormat','struct');
longitude = cellfun(@(m) double(m.Longitude), msgs);
latitude = cellfun(@(m) double(m.Latitude), msgs);
altitude = cellfun(@(m) double(m.Altitude), msgs);

l0 = longitude(1);
mu0 = latitude(1);
h_ref = altitude(1);


% for i = 1:length(altitude)
%     l = longitude(i);
%     mu = latitude(i);
%     h = altitude(i);
%     


% figure
% plot(ts, 'LineWidth', 1)
