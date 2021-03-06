function [p_n_nb_data,v_n_nb_data,q_nb_data, bacc_b_nb_data, bars_b_nb_data, f_b_imu_data, omega_b_imu_data,time, g_n_nb, v_abs, acc_std, ars_std] = SkidPadSim(simtime ,frequency, enable_plots)
%CIRCLESIM Summary of this function goes here
%   Radius = 50 [m]
%   Absolute velocity = (2 * pi * R)/30 [m/s] (30 sec per circle)

rad2deg = 180/pi;
deg2rad = pi/180;

h = 1/frequency;
N = simtime/h;

t_goal = 5.5;

% Circle data
R = 9.125;
v_abs = (2 * pi * R) / t_goal;
% inclination = 20;
% height = 18.20;

g = 9.81;
t = 0;
count = 1;


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


% Initialize
p_n_nb = [0 0 0]';
v_n_nb = [0 0 0]';
% att_n_nb = [0 0 0]';

roll = deg2rad * 10;
pitch = deg2rad * 0;
yaw = deg2rad * 0;
q_nb = euler2q(roll, pitch, yaw);
q_nb = q_nb / norm(q_nb);
% q_nb = [1 0 0 0]';

% Bias
acc_bias = [-0.4 -0.5 0.3]';
ars_bias = [-.030 0.02 -.02]';
acc_bias_std = (0.04 * 0.001 * 9.80665) * sqrt(h);
ars_bias_std = ((10 * deg2rad) / 3600) * sqrt(h);
bacc_b_nb = acc_bias;
bars_b_nb = ars_bias;

% Standard deviations (for computation of noise)
acc_std = (0.14 * 0.001 * 9.80665) / sqrt(h);
ars_std = (0.0035 * deg2rad) / sqrt(h);

v_b_nb = [0 0 0]';

% alpha = 0;
% pitch_rate = deg2rad * 11.31*cos(alpha - (pi/2));
% roll_rate = deg2rad * 11.31*cos(alpha);
% yaw_rate = pi/15;
% omega_b_nb = [pitch_rate roll_rate yaw_rate]';

% omega_b_nb = [0 0 pi/(5.5/2)]';
% a_b_nb = [0 (v_abs^2)/R 0]';

omega_b_nb = [0 0 0]';
a_b_nb = [0 0 0]';

g_n_nb = [0 0 g]';

skid_flag = false;
start = false;

for k = 1:N
    
    time(k) = t;
    
    if (t > 14*t_goal) && (start == false)
        omega_b_nb = [0 0 0]';
        a_b_nb = [v_abs/t_goal 0 0]';
        start = true;
    end
    
    
    if (t > 15*t_goal) && (skid_flag == false)
        omega_b_nb = [0 0 pi/(t_goal/2)]';
        a_b_nb = [0 (v_abs^2)/R 0]';
        skid_flag = true;
    end
    
    modulo = mod(round(t,5),(t_goal/2)); 
    if (modulo == 0)&& (skid_flag == true)
        count = count + 1;
        if ( count == 4 )
            omega_b_nb = -omega_b_nb;
            a_b_nb = -a_b_nb;
            count = 0; 
        end
    end
    p_n_nb_data(:,k) = p_n_nb;
    v_n_nb_data(:,k) = v_n_nb;
%     att_n_nb_data(:,k) = att_n_nb;
    q_nb_data(:,k) = q_nb;
    bacc_b_nb_data(:,k) = acc_bias;
    bars_b_nb_data(:,k) = ars_bias;
    
%     phi = att_n_nb_data(1,k);
%     theta = att_n_nb_data(2,k);
%     psi = att_n_nb_data(3,k);

%     [J,R_nb,T_nb] = eulerang(phi, theta, psi);

    [J,R_nb, T_nb] = quatern(q_nb_data(:,k));
%     rng(working_seed);
    acc_noise = acc_std * randn(3,1);
    ars_noise = ars_std * randn(3,1);
%     acc_noise = 0.001 * wgn(3, 1, 1);
%     ars_noise = 0.0005 * wgn(3, 1, 1); 

    f_b_imu = a_b_nb + Smtrx(omega_b_nb)*v_b_nb - (R_nb')*g_n_nb + bacc_b_nb + acc_noise; 
    omega_b_imu = omega_b_nb + bars_b_nb + ars_noise;
    q_omega = qbuild((omega_b_imu - bars_b_nb - ars_noise),h);
    
    f_b_imu_data(:,k) = f_b_imu;
    omega_b_imu_data(:,k) = omega_b_imu;
       
    a_b_nb = (f_b_imu - bacc_b_nb - acc_noise) + (R_nb') * g_n_nb;
    p_n_nb = p_n_nb + (h * v_n_nb) + (0.5 * h * h * R_nb * a_b_nb);
    v_n_nb = v_n_nb + (h * R_nb * a_b_nb);
%     att_n_nb = att_n_nb + (h * omega_b_imu);
%     att_n_nb = ssa(att_n_nb,'rad'); 
    q_nb = quatprod(q_nb,q_omega);
    q_nb = q_nb / norm(q_nb);
    bacc_b_nb = acc_bias + acc_bias_std * randn(3, 1);
    bars_b_nb = ars_bias + ars_bias_std * randn(3, 1);

    t = t + h;
%     t_count = t_count + h;
    
%     alpha = alpha + (pi/(15*h));
%     [ax_b_nb, ay_b_nb, az_b_nb] = trueBodyAccel(alpha, inclination);
%     a_b_nb = [ax_b_nb, ay_b_nb, az_b_nb]';
%     [omegax_b_nb, omegay_b_nb, omegaz_b_nb] = trueBodyAngRate(alpha, height);
%     omega_b_nb = [omegax_b_nb, omegay_b_nb, omegaz_b_nb]';
    
    
%     if(t_count >= 15)
%        t_count = 0;
%        omega_b_nb = [ -omega_b_nb(1) ; -omega_b_nb(2) ; omega_b_nb(3) ];
%     end

end


% FUNCTIONS
function [ax_b_nb, ay_b_nb, az_b_nb] = trueBodyAccel(alpha, height)
       ax_b_nb = 0;
       ay_b_nb = (v_abs^2)/R;
       az_b_nb = (height/2)*sin(alpha-(pi/2)); 
end

function [omegax_b_nb, omegay_b_nb, omegaz_b_nb] = trueBodyAngRate(alpha, inclination)
       omegax_b_nb = (pi/180) * inclination*cos(alpha - (pi/2));
       omegay_b_nb = (pi/180) * inclination*cos(alpha);
       omegaz_b_nb = pi/15;
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

%     % Attitude
%     figure(3)
%     figure(gcf);
%     subplot(3, 1, 1)
%     hold on;
%     plot(time, rad2deg*q2euler(q_nb_data(1,:)), 'Color', 'black', 'Linewidth', 1.5);
%     ylabel('Roll [m]')
%     legend('True');
%     title('Attitude');
% 
%     subplot(3, 1, 2)
%     hold on;
%     plot(time, rad2deg*q2euler(q_nb_data(2,:)), 'Color', 'black', 'Linewidth', 1.5);
%     ylabel('Pitch [m]')
%     legend('True');
% 
%     subplot(3, 1, 3)
%     hold on;
%     plot(time, rad2deg*q2euler(q_nb_data(3,:)), 'Color', 'black', 'Linewidth', 1.5);
%     xlabel('Time [s]');
%     ylabel('Yaw [m]')
%     legend('True');

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
