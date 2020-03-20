clear all;
clc;
close all;

deg2rad = pi/180;   
rad2deg = 180/pi;

Z3 = zeros(3,3);
I3 = eye(3);

simtime = 300;
f_samp  = 100;          %imu frequency
f_low   = 10;           %aiding frequency
h       = 1/f_samp;     %sampling time
N       = simtime/h;    %number of iterations

%init values
delta_x = zeros(15, 1);    % [dp, dv, dbacc, dtheta, dbars]
x_ins = zeros(16, 1);      % [p, v, bacc, q, bars]
y = zeros(6,1);

std_pos = 2;
std_att = 5 * deg2rad;

phi_ins = 0;
theta_ins = 0;
psi_ins = 0;

q_n_ins = euler2q(phi_ins,theta_ins,psi_ins);
q_n_ins = q_n_ins/norm(q_n_ins);


% data storage
time_data = zeros(1, N);
ins_data = zeros(15, N);
att_n_nb = zeros(3, N);
g_err_data = zeros(3,N);

% from sim
[seed, p_n_nb, v_n_nb, q_nb, bacc_b_nb, bars_b_nb, f_b_imu, omega_b_imu,time, g_n_nb, v_abs] = CircleSim(simtime,f_samp,0);
% [p_n_nb, v_n_nb, att_n_nb, f_b_imu, w_b_imu,time] = StandStillSim(simtime,f_samp,0);

% dual gnss config
r_b_1 = [0.5 0 0]';
r_b_2 = [-0.5 0 0.3]';

%initialization of kalman filter
f_b_imu_0 = [0 0 0]';
omega_b_imu_0 = [0 0 0]';
bacc_b_ins = [0 0 0]';
bars_b_ins = [0 0 0]';
E_prev = zeros(15,12);
ErrorStateKalman_sola(r_b_1, r_b_2, E_prev,0, 0, f_low, 1, f_b_imu_0, omega_b_imu_0, g_n_nb, bacc_b_ins, bars_b_ins);

% init
x_ins(1:3) = [5;6;7]; % for testing av ESKF
x_ins(4:6) = [0;0;0]; % for testing av ESKF
x_ins(10:13) = [1 0 0 0]'; % for testing av ESKF
% x_ins(4:6) = [v_abs;0;0]; % korrekt initialisering. For testing av INS

% matrices
C_ins = [I3 Z3 Z3 Z3 Z3
         Z3 Z3 Z3 I3 Z3];
     
     
count = 10;

for k = 1:N
   t = k * h;
   time_data(k) = t;
   g_err_data(1:3,k) = g_n_nb;
   
   % split state vector
   p_n_ins = x_ins(1:3);
   v_n_ins = x_ins(4:6);
   bacc_b_ins = x_ins(7:9);
   q_n_ins = x_ins(10:13);
   q_n_ins = q_n_ins/norm(q_n_ins);
   bars_b_ins = x_ins(14:16);
   
   % compute rotation vector
   R_nb_ins = Rquat(q_n_ins);
   
   % store current state 
   [phi_ins, theta_ins, psi_ins] = q2euler(q_n_ins);
   att_n_ins = [phi_ins, theta_ins, psi_ins]';
   ins_data(:,k) = [p_n_ins ; v_n_ins ; bacc_b_ins ; att_n_ins ; bars_b_ins]; 
   
   [phi_t,theta_t,psi_t] = q2euler(q_nb(:,k));
   att_n_nb(1:3,k) = [phi_t theta_t psi_t]';
   
   % compute acceleration and angular rate
   a_n_ins = R_nb_ins*(f_b_imu(:,k) - bacc_b_ins) + g_n_nb;
   omega_b_ins = omega_b_imu(:,k) - bars_b_ins;
   
   % compute quaternion from angular rates 
   q_omega_b_ins = qbuild(omega_b_ins, h);
   
   % update nominal states with imu input
   p_n_ins = p_n_ins + (h * v_n_ins) + (0.5 * h * h * a_n_ins);
   v_n_ins = v_n_ins +  (h * a_n_ins);
%    bacc_b_ins = bacc_b_ins;
   q_n_ins = quatprod(q_n_ins, q_omega_b_ins); 
   q_n_ins = q_n_ins/norm(q_n_ins);
%    bars_b_ins = bars_b_ins;
   
   x_ins = [p_n_ins ; v_n_ins ; bacc_b_ins ; q_n_ins ; bars_b_ins];
   
   count = count + 1;

   if (count >= 10)
        count = 0;
        
        % noisy measurements
        p_meas = p_n_nb(1:3,k) +  0.001 * wgn(3, 1, 1);
        q_meas = q_nb(1:4,k) + 0.00005 * wgn(4, 1, 1);
        
        q_conj = quatconj(q_n_ins')';
        delta_q = quatprod(q_conj, q_meas);
        delta_theta = 2*delta_q(2:4);
        
        p_gnss_1 = p_n_nb(1:3,k) + I3*R_nb_ins*r_b_1 - Smtrx(R_nb_ins*r_b_1)*delta_theta;
        p_gnss_2 = p_n_nb(1:3,k) + I3*R_nb_ins*r_b_2 - Smtrx(R_nb_ins*r_b_2)*delta_theta;
        
        R_nb_t = Rquat(q_nb(1:4,k)');
        
        p_gnss_1 = p_n_nb(1:3,k) + R_nb_t * r_b_1;
        p_gnss_2 = p_n_nb(1:3,k) + R_nb_t * r_b_2;
%         p_gnss_3 = p_n_nb(1:3,k) + R_nb_t * r_b_3;
        
        
        p_hat_1 = p_n_ins + R_nb_ins * r_b_1;
        p_hat_2 = p_n_ins + R_nb_ins * r_b_2;
%         p_hat_3 = p_n_ins + R_nb_ins * r_b_3;
        
        
%         delta_y = [(p_gnss_1 - p_hat_1) ; (p_gnss_2 - p_hat_2); (p_gnss_3 - p_hat_3); delta_theta]; 
        delta_y = [(p_gnss_1 - p_hat_1) ; (p_gnss_2 - p_hat_2)];
        % compute error state with ESKF
        [delta_x, E_prev] = ErrorStateKalman_sola(r_b_1, r_b_2, E_prev, delta_y, R_nb_ins, f_low, 0, f_b_imu, omega_b_imu, g_n_nb, bacc_b_ins, bars_b_ins);

% 
%         % delta y
%         delta_p = p_meas - p_n_ins; 
% 
%         q_conj = quatconj(q_n_ins')';
%         delta_q = quatprod(q_conj, q_meas);
%         delta_theta = 2*delta_q(2:4);
%         
%         delta_y = [delta_p; delta_theta];         
%         
%         % compute error state with ESKF
%         [delta_x, E_prev] = ErrorStateKalman_sola(E_prev, delta_y, R_nb_ins, f_low, 0, f_b_imu, omega_b_imu, g_n_nb, bacc_b_ins, bars_b_ins);
%         
        % inject error into nominal state
        x_ins(1:9) = x_ins(1:9) + delta_x(1:9);
        x_ins(14:16) = x_ins(14:16) + delta_x(13:15);
        
%         g_n_nb = g_n_nb + delta_x(16:18);
        
        h_low = 1/10;
        q_delta_omega = qbuild(delta_x(10:12)/h_low, h_low);
        x_ins(10:13) = quatprod(x_ins(10:13), q_delta_omega);
        x_ins(10:13) = x_ins(10:13)/norm(x_ins(10:13)); 
   end
    
   
end


% PLOTS
   
% POSITION
figure(1)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, ins_data(1,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, p_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
ylabel('X position [m]')
legend('Est', 'True');
title('Position');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, ins_data(2,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, p_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Y position [m]')
legend('Est', 'True');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, ins_data(3,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, p_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Z position [m]')
legend('Est', 'True');
grid on;
saveas(gcf,'results/Position.jpeg')

% POSITION ERROR
figure(2)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, p_n_nb(1,:) - ins_data(1,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('X position [m]')
legend('Est');
title('Position Error');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, p_n_nb(2,:) - ins_data(2,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('Y position [m]')
legend('Est');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, p_n_nb(3,:) - ins_data(3,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('Z position [m]')
legend('Est');
grid on;
saveas(gcf,'results/PositionError.jpeg')


% VELOCITIES
figure(3)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, ins_data(4,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, v_n_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('X velocity [m]')
legend('Est', 'True');
title('Velocity');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, ins_data(5,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, v_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Y velocity [m]')
legend('Est', 'True');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, ins_data(6,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, v_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Z velocity [m]')
legend('Est', 'True');
grid on;
saveas(gcf,'results/Velocity.jpeg')

% VELOCITY ERROR
figure(4)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, v_n_nb(1,:) - ins_data(4,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('X velocity error [m/s]')
legend('Est');
title('Velocity Error');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, v_n_nb(2,:) - ins_data(5,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('Y velocity error [m/s]')
legend('Est');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, v_n_nb(3,:) - ins_data(6,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('Z velocity error [m/s]')
legend('Est');
grid on;
saveas(gcf,'results/VelocityError.jpeg')



% ACCEL BIAS
figure(5)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, ins_data(7,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, bacc_b_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('X acc bias [deg]')
legend('Est', 'True');
title('Accelerometer bias');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, ins_data(8,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, bacc_b_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Y acc bias [deg]')
legend('Est', 'True');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, ins_data(9,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, bacc_b_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Z acc bias [deg]')
legend('Est', 'True');
grid on;
saveas(gcf,'results/AccelBiast.jpeg')


% ACCEL BIAS ERROR
figure(6)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, bacc_b_nb(1,:) - ins_data(7,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('X accel bias error [m/s^2]')
legend('Est');
title('Accelerometer Bias Error');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, bacc_b_nb(2,:) - ins_data(8,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('Y accel bias error [m/s^2]')
legend('Est');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, bacc_b_nb(3,:) - ins_data(9,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('Z accel bias error [m/s^2]')
legend('Est');
grid on;
saveas(gcf,'results/AccelBiasError.jpeg')



% ATTITUDE 
figure(7)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, rad2deg*ins_data(10,:),'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, rad2deg*att_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
ylabel('Roll angle [deg]')
legend('Est', 'True');
title('Attitude');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, rad2deg*ins_data(11,:), 'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, rad2deg*att_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Pitch angle [deg]')
legend('Est', 'True');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, rad2deg*ins_data(12,:), 'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, rad2deg*att_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
xlabel('Time [s]');
ylabel('yaw angle [deg]')
legend('Est', 'True');
grid on;
saveas(gcf,'results/Attitude.jpeg')


% ATTITUDE ERROR
figure(8)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, ssa(rad2deg*att_n_nb(1,:), 'deg') - ssa(rad2deg*ins_data(10,:),'deg') ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('X velocity error [m/s]')
legend('Est');
title('Attotide Error');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, ssa(rad2deg*att_n_nb(2,:),'deg') - ssa(rad2deg*ins_data(11,:), 'deg') ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('Y velocity error [m/s]')
legend('Est');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, ssa(rad2deg*att_n_nb(3,:),'deg') - ssa(rad2deg*ins_data(12,:),'deg') ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('Z velocity error [m/s]')
legend('Est');
grid on;
saveas(gcf,'results/AttitudeError.jpeg')



% GYRO BIAS
figure(9)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, ins_data(13,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, bars_b_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Roll bias [deg]')
legend('Est', 'True');
title('Gyro bias');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, ins_data(14,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, bars_b_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Pitch bias [deg]')
legend('Est', 'True');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, ins_data(15,:) ,'Color', [1,165/255, 0], 'Linewidth', 2);
plot(time, bars_b_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Yaw bias [deg]')
legend('Est', 'True');
grid on;
saveas(gcf,'results/GyroBias.jpeg')

% GYRA BIAS ERROR
figure(10)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, bars_b_nb(1,:) - ins_data(13,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('X gyro bias error [rad/s^2]')
legend('Est');
title('Gyro Bias Error');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, bars_b_nb(2,:) - ins_data(14,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('Y gyro bias error [rad/s^2]')
legend('Est');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, bars_b_nb(3,:) - ins_data(15,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylabel('Z gyro bias error [rad/s^2]')
legend('Est');
grid on;
saveas(gcf,'results/GyroBiasError.jpeg')


% % Accelerometer input
% figure(11)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, f_b_imu(1,:),'Color', 'black', 'Linewidth', 1.5);
% ylabel('f_b_imu X [m/s^2]')
% legend('Measured');
% title('Accelerometer input');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, f_b_imu(2,:),'Color', 'black', 'Linewidth', 1.5);
% ylabel('f_b_imu y [m/s^2]')
% legend('Measured');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, f_b_imu(3,:),'Color', 'black', 'Linewidth', 1.5);
% ylabel('f_b_imu Z [m/s^2]')
% legend('Measured');
% grid on;
% % saveas(gcf,'AccelInput.jpeg')
% 
% 
% 
% % Gyrometer input
% figure(12)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, rad2deg*omega_b_imu(1,:),'Color', 'black', 'Linewidth', 1.5);
% ylabel('omega_b_imu X [m/s^2]')
% legend('Measured');
% title('Gyrometer input');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, rad2deg*omega_b_imu(2,:),'Color', 'black', 'Linewidth', 1.5);
% ylabel('omega_b_imu y [m/s^2]')
% legend('Measured');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, rad2deg*omega_b_imu(3,:),'Color', 'black', 'Linewidth', 1.5);
% ylabel('omega_b_imu Z [m/s^2]')
% legend('Measured');
% grid on;
% saveas(gcf,'GyroInput.jpeg')



% gravity error
figure(13)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, g_err_data(1,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('gravity err X [m/s^2]')
legend('Estimated');
title('Gravity Error');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, g_err_data(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('gravity err Y [m/s^2]')
legend('Estimated');
grid on;


subplot(3, 1, 3)
hold on;
plot(time, g_err_data(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('gravity err Z [m/s^2]')
legend('Estimated');
grid on;

saveas(gcf,'results/GravityErrorState.jpeg')

% Position map
figure(14)
figure(gcf)
subplot(1, 1, 1)
plot(ins_data(2,:), ins_data(1,:), 'Color', 'black', 'Linewidth', 1.5);
xlabel('Y position [m]');
ylabel('X position [m]');
title('Position Plot');
grid on;
saveas(gcf,'results/PositionMap.jpeg')

