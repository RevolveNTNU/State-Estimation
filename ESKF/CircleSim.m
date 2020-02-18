function [p_n_nb,v_n_nb,att_n_nb] = CircleSim(simtime ,frequency)
%CIRCLESIM Summary of this function goes here
%   Radius = 50 [m]
%   Absolute velocity = (2 * pi * R)/30 [m/s] (30 sec per circle)

rad2deg = 180/pi;

h = 1/frequency;
N = simtime/h;

% Circle data
R = 50;
v_abs = (2 * pi * R) / 30; 

g = 9.81;
t = 0;

% Allocate data
time = zeros(1, N);
p_n_nb = zeros(3, N);
v_n_nb = zeros(3, N);
att_n_nb = zeros(3, N);
acc_n_nb = zeros(3, N);
angvel_b_nb = zeros(3, N);

% Initialize
p_n_nb(1:3,1) = [0 0 0]';
v_n_nb(1:3,1) = [v_abs 0 0]';
att_n_nb(1:3,1) = [0 0 0]';


v_b_nb = [0 0 0]';
w_b_nb = [0 0 pi/15]';
a_b_nb = [0 (v_abs^2)/R 0]';


for i = 2:N

    phi = att_n_nb(1,i-1);
    theta = att_n_nb(2,i-1);
    psi = att_n_nb(3,i-1);

    [J,R_nb,T_nb] = eulerang(phi, theta, psi);

    f_b_imu = a_b_nb + Smtrx(w_b_nb)*v_b_nb; % - (R_nb')*[0 0 g]'; 
    w_b_imu = w_b_nb;

    p_n_nb(1:3, i) = p_n_nb(1:3, i-1) + (h * v_n_nb(1:3, i-1)) + (0.5 * h * h * f_b_imu);
    v_n_nb(1:3, i) = v_n_nb(1:3, i-1) + (h * R_nb * f_b_imu);
    att_n_nb(1:3, i) = att_n_nb(1:3, i-1) + (h * w_b_imu);

    time(i) = t;

    t = t + h;


end

% PLOTS


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


