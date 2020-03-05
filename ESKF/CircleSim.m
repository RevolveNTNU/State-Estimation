function [p_n_nb_data,v_n_nb_data,q_nb_data, bacc_b_nb_data, bars_b_nb_data, f_b_imu_data, omega_b_imu_data,time, g_n_nb, v_abs] = CircleSim(simtime ,frequency, enable_plots)
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
p_n_nb_data = zeros(3, N);
v_n_nb_data = zeros(3, N);
% att_n_nb_data = zeros(3, N);
q_nb_data = zeros(4,N); 
bacc_b_nb_data = zeros(3, N);
bars_b_nb_data = zeros(3, N);
f_b_imu_data = zeros(3,N);
omega_b_imu_data = zeros(3,N);

% acc_n_nb = zeros(3, N);
% angvel_b_nb = zeros(3, N);

% Initialize
p_n_nb = [0 0 0]';
v_n_nb = [v_abs 0 0]';
% att_n_nb = [0 0 0]';
q_nb = [1 0 0 0]';

acc_bias = 0.001*[-0.4 -0.5 0.3]';
ars_bias = 0.001*[-.030 0.02 -.02]';
bacc_b_nb = acc_bias;
bars_b_nb = ars_bias;


v_b_nb = [0 0 0]';
omega_b_nb = [0 0 pi/15]';
a_b_nb = [0 (v_abs^2)/R 0]';
g_n_nb = [0 0 g]';


for k = 1:N
    
    time(k) = t;
    p_n_nb_data(:,k) = p_n_nb;
    v_n_nb_data(:,k) = v_n_nb;
%     att_n_nb_data(:,k) = att_n_nb;
    q_nb_data(:,k) = q_nb;
    bacc_b_nb_data(:,k) = bacc_b_nb;
    bars_b_nb_data(:,k) = bars_b_nb;
    
%     phi = att_n_nb_data(1,k);
%     theta = att_n_nb_data(2,k);
%     psi = att_n_nb_data(3,k);

%     [J,R_nb,T_nb] = eulerang(phi, theta, psi);

    [J,R_nb, T_nb] = quatern(q_nb_data(:,k));

    f_b_imu = a_b_nb + Smtrx(omega_b_nb)*v_b_nb - (R_nb')*g_n_nb + bacc_b_nb; 
    omega_b_imu = omega_b_nb + bars_b_nb;
    q_omega = qbuild((omega_b_imu - bars_b_nb),h);
    
    f_b_imu_data(:,k) = f_b_imu;
    omega_b_imu_data(:,k) = omega_b_imu;
       
    a_b_nb = f_b_imu - bacc_b_nb + (R_nb') * g_n_nb;
    p_n_nb = p_n_nb + (h * v_n_nb) + (0.5 * h * h * R_nb * a_b_nb);
    v_n_nb = v_n_nb + (h * R_nb * a_b_nb);
%     att_n_nb = att_n_nb + (h * omega_b_imu);
%     att_n_nb = ssa(att_n_nb,'rad'); 
    q_nb = quatprod(q_nb,q_omega);
    bacc_b_nb = acc_bias + 0.001 * wgn(3, 1, 1);
    bars_b_nb = ars_bias + 0.00005 * wgn(3, 1, 1);

    t = t + h;

end

% PLOTS

if (enable_plots)
    % Position
    figure(1)
    figure(gcf);
    subplot(3, 1, 1)
    hold on;
    plot(time, p_n_nb_data(1,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('X position [m]')
    legend('True');
    title('Position');

    subplot(3, 1, 2)
    hold on;
    plot(time, p_n_nb_data(2,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('Y position [m]')
    legend('True');

    subplot(3, 1, 3)
    hold on;
    plot(time, p_n_nb_data(3,:), 'Color', 'black', 'Linewidth', 1.5);
    xlabel('Time [s]');
    ylabel('Z position [m]')
    legend('True');


    % Velocity
    figure(2)
    figure(gcf);
    subplot(3, 1, 1)
    hold on;
    plot(time, v_n_nb_data(1,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('X velocity [m]')
    legend('True');
    title('Velocity');

    subplot(3, 1, 2)
    hold on;
    plot(time, v_n_nb_data(2,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('Y velocity [m]')
    legend('True');

    subplot(3, 1, 3)
    hold on;
    plot(time, v_n_nb_data(3,:), 'Color', 'black', 'Linewidth', 1.5);
    xlabel('Time [s]');
    ylabel('Z velocity [m]')
    legend('True');

    % Attitude
    figure(3)
    figure(gcf);
    subplot(3, 1, 1)
    hold on;
    plot(time, rad2deg*att_n_nb_data(1,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('Roll [m]')
    legend('True');
    title('Attitude');

    subplot(3, 1, 2)
    hold on;
    plot(time, rad2deg*att_n_nb_data(2,:), 'Color', 'black', 'Linewidth', 1.5);
    ylabel('Pitch [m]')
    legend('True');

    subplot(3, 1, 3)
    hold on;
    plot(time, rad2deg*att_n_nb_data(3,:), 'Color', 'black', 'Linewidth', 1.5);
    xlabel('Time [s]');
    ylabel('Yaw [m]')
    legend('True');

    % Position map
    figure(6)
    figure(gcf)
    subplot(1, 1, 1)
    plot(-p_n_nb_data(2,:), p_n_nb_data(1,:), 'Color', 'black', 'Linewidth', 1.5);
    xlabel('Y position [m]');
    ylabel('X position [m]');
    title('Position Plot');
end
    
end


