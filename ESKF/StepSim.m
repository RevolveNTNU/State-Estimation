function [p_n_nb,v_n_nb,att_n_nb,f_b_imu, w_b_imu,time] = StepSim(simtime ,frequency, enable_plots)
%STANDSTILLSIM Summary of this function goes here
%   Detailed explanation goes here

deg2rad = pi/180;
rad2deg = 180/pi;

h = 1/frequency;
N = simtime/h;

t = 0;

% Allocate data
time = zeros(1, N);
p_n_nb = zeros(3, N);
v_n_nb = zeros(3, N);
att_n_nb = zeros(3, N);

for i = 2:N
    t = t + h;
   time(i) = t;
end

const_pos_1 = [100,100,0]';
const_pos_2 = [200,150,0]';
const_vel = [0 0 0]';
const_att = [0,0,deg2rad*50]';

p_n_nb = repmat(const_pos_1, 1, N);
for i = N/2:N
    p_n_nb(1:3, i) = const_pos_2;
end
v_n_nb = repmat(const_vel, 1, N);
att_n_nb = repmat(const_att, 1, N);
f_b_imu = [0 0 0]';
w_b_imu = [0 0 0]';


if (enable_plots)
    % Position
    figure(1)
    figure(gcf);
    subplot(3, 1, 1)
    hold on;
    plot(time, p_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('X position [m]')
    legend('True');
    title('Position');

    subplot(3, 1, 2)
    hold on;
    plot(time, p_n_nb(2,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('Y position [m]')
    legend('True');

    subplot(3, 1, 3)
    hold on;
    plot(time, p_n_nb(3,:), 'Color', 'black', 'Linewidth', 1.5);
    xlabel('Time [s]');
    ylabel('Z position [m]')
    legend('True');


    % Velocity
    figure(2)
    figure(gcf);
    subplot(3, 1, 1)
    hold on;
    plot(time, v_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('X velocity [m]')
    legend('True');
    title('Velocity');

    subplot(3, 1, 2)
    hold on;
    plot(time, v_n_nb(2,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('Y velocity [m]')
    legend('True');

    subplot(3, 1, 3)
    hold on;
    plot(time, v_n_nb(3,:), 'Color', 'black', 'Linewidth', 1.5);
    xlabel('Time [s]');
    ylabel('Z velocity [m]')
    legend('True');

    % Attitude
    figure(3)
    figure(gcf);
    subplot(3, 1, 1)
    hold on;
    plot(time, rad2deg*att_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('Roll [m]')
    legend('True');
    title('Attitude');

    subplot(3, 1, 2)
    hold on;
    plot(time, rad2deg*att_n_nb(2,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('Pitch [m]')
    legend('True');

    subplot(3, 1, 3)
    hold on;
    plot(time, rad2deg*att_n_nb(3,:), 'Color', 'black', 'Linewidth', 1.5);
    xlabel('Time [s]');
    ylabel('Yaw [m]')
    legend('True');

    % Position map
    figure(6)
    figure(gcf)
    subplot(1, 1, 1)
    plot(-p_n_nb(2,:), p_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
    xlabel('Y position [m]');
    ylabel('X position [m]');
    title('Position Plot');
end
 
end

