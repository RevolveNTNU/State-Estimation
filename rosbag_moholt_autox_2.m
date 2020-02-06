clear all;
clc;
close all;

bag = rosbag('C:\Users\tonja\OneDrive\Documents\rosbags\moholt_autox_2.bag');
bag.AvailableTopics;
bag.MessageList;
start = bag.StartTime;
bagselect = select(bag,'Topic','/vn300/Imu');

% ts = timeseries(bagselect, 'LinearAcceleration.X', 'LinearAcceleration.Y', 'LinearAcceleration.Z' );
% ts.Data;

% msgs = readMessages(bagselect);
msgs = readMessages(bagselect, 1:4:12000);
size(msgs);
msgs{3};
ax = cellfun(@(m) double(m.LinearAcceleration.X), msgs);
ay = cellfun(@(m) double(m.LinearAcceleration.Y), msgs);
az = cellfun(@(m) double(m.LinearAcceleration.Z), msgs);

omegax = cellfun(@(m) double(m.AngularVelocity.X), msgs);
omegay = cellfun(@(m) double(m.AngularVelocity.Y), msgs);
omegaz = cellfun(@(m) double(m.AngularVelocity.Z), msgs);

% figure
% plot(ts, 'LineWidth', 1)
