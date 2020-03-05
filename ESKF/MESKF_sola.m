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

q_b_ins = euler2q(phi_ins,theta_ins,psi_ins);
q_b_ins = q_b_ins/norm(q_b_ins);


% data storage
time_data = zeros(1, N);
ins_data = zeros(15, N);
att_n_nb = zeros(3, N);

% from sim
[p_n_nb, v_n_nb, q_nb, bacc_b_nb, bars_b_nb, f_b_imu, omega_b_imu,time, g_n_nb, v_abs] = CircleSim(simtime,f_samp,0);
% [p_n_nb, v_n_nb, att_n_nb, f_b_imu, w_b_imu,time] = StandStillSim(simtime,f_samp,0);

%initialization of kalman filter
f_b_imu_0 = [0 0 0]';
omega_b_imu_0 = [0 0 0]';
bacc_b_ins = [0 0 0]';
bars_b_ins = [0 0 0]';
ErrorStateKalman_sola(0, 0, f_low, 1, f_b_imu_0, omega_b_imu_0, g_n_nb, bacc_b_ins, bars_b_ins);

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
   
   % split state vector
   p_n_ins = x_ins(1:3);
   v_n_ins = x_ins(4:6);
   bacc_b_ins = x_ins(7:9);
   q_b_ins = x_ins(10:13);
   q_b_ins = q_b_ins/norm(q_b_ins);
   bars_b_ins = x_ins(14:16);
   
   % compute rotation vector
   R_nb_ins = Rquat(q_b_ins);
   
   % store current state 
   [phi_ins, theta_ins, psi_ins] = q2euler(q_b_ins);
   att_n_ins = [phi_ins, theta_ins, psi_ins]';
   ins_data(:,k) = [p_n_ins ; v_n_ins ; bacc_b_ins ; att_n_ins ; bars_b_ins]; 
   
   [phi_t,theta_t,psi_t] = q2euler(q_nb(:,k));
   att_n_nb(1:3,k) = [phi_t theta_t psi_t]';
   
   % compute acceleration and angular rate
   a_n_ins = R_nb_ins*(f_b_imu(:,k) - bacc_b_ins) + g_n_nb;
   omega_b_ins = omega_b_imu(:,k) - bars_b_ins;
   
   % compute quaternion from angular rates 
   q_omega_b_ins = qbuild(omega_b_ins, h);
%    alpha1 = omega_b_ins(1)*h;
%    alpha2 = omega_b_ins(2)*h;
%    alpha3 = omega_b_ins(3)*h;
%    
%    q_omega = [cos(alpha1/2) ; sin(alpha1/2)*[1 0 0]'] + [cos(alpha2/2) ; sin(alpha2/2)*[0 1 0]'] + [cos(alpha3/2) ; sin(alpha3/2)*[0 0 1]'];  
%    q_omega = quatprod( [cos(alpha1/2) ; sin(alpha1/2)*[1 0 0]'], [cos(alpha2/2) ; sin(alpha2/2)*[0 1 0]']);
%    q_omega = quatprod( q_omega, [cos(alpha3/2) ; sin(alpha3/2)*[0 0 1]']);
%    q_omega = q_omega/norm(q_omega);

   
   % update nominal states with imu input
   p_n_ins = p_n_ins + (h * v_n_ins) + (0.5 * h * h * a_n_ins);
   v_n_ins = v_n_ins +  (h * a_n_ins);
   bacc_b_ins = bacc_b_ins;
   q_b_ins = quatprod(q_b_ins, q_omega_b_ins); %quatprod(q_ins, [0 ; omega_b_ins]);
   q_b_ins = q_b_ins/norm(q_b_ins);
   bars_b_ins = bars_b_ins;
   
   x_ins = [p_n_ins ; v_n_ins ; bacc_b_ins ; q_b_ins ; bars_b_ins];
   
   count = count + 1;

   if (count >= 10)
        count = 0;
        
        % perfect measurement
        y(1:3) = p_n_nb(1:3,k); % + (std_pos * randn(1) * ones(1, 3))';
        y(4:6) = att_n_nb(1:3,k); % + (std_pos * randn(1) * ones(1, 3))';
        
%         y_ins = [x_ins(1:3) ; x_ins(11:13)]
        [phi_ins, theta_ins, psi_ins] = q2euler(q_b_ins);
        y_ins = [ p_n_ins ; [phi_ins, theta_ins, psi_ins]'];
        
        delta_y = y - y_ins;
        
        % compute error state with ESKF
        delta_x = ErrorStateKalman_sola(delta_y, R_nb_ins, f_low, 0, f_b_imu, omega_b_imu, g_n_nb, bacc_b_ins, bars_b_ins);
        
        % add error to nominal state
        x_ins(1:9) = x_ins(1:9) + delta_x(1:9);
        x_ins(14:16) = x_ins(14:16) + delta_x(13:15);
        
        % multiplicative part
        h_low = 1/10;
        q_delta_omega = qbuild(delta_x(10:12)/h_low, h_low);
        
%         alpha1 = delta_x(10);
%         alpha2 = delta_x(11);
%         alpha3 = delta_x(12);
%         q_omega = [cos(alpha1/2) ; sin(alpha1/2)*[1 0 0]'] + [cos(alpha2/2) ; sin(alpha2/2)*[0 1 0]'] + [cos(alpha3/2) ; sin(alpha3/2)*[0 0 1]'];  
%         q_omega = quatprod( [cos(alpha1/2) ; sin(alpha1/2)*[1 0 0]'], [cos(alpha2/2) ; sin(alpha2/2)*[0 1 0]']);
%         q_omega = quatprod( q_omega, [cos(alpha3/2) ; sin(alpha3/2)*[0 0 1]']);
%         q_omega = q_omega/norm(q_omega);
   
%         delta_q = euler2q(delta_x(10), delta_x(11), delta_x(12));
%         delta_q = [1 ; 0.5 * delta_x(10:12)];
%         delta_q = delta_q / norm(delta_q);
        
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
plot(time, ins_data(1,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, p_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
ylabel('X position [m]')
legend('Est', 'True');
title('Position');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, ins_data(2,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, p_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Y position [m]')
legend('Est', 'True');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, ins_data(3,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, p_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Z position [m]')
legend('Est', 'True');
grid on;

% VELOCITIES
figure(2)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, ins_data(4,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, v_n_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('X velocity [m]')
legend('Est', 'True');
title('Velocity');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, ins_data(5,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, v_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Y velocity [m]')
legend('Est', 'True');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, ins_data(6,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, v_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Z velocity [m]')
legend('Est', 'True');
grid on;

% ACCEL BIAS
figure(3)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, ins_data(7,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, bacc_b_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('X acc bias [deg]')
legend('Est', 'True');
title('Accelerometer bias');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, ins_data(8,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, bacc_b_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Y acc bias [deg]')
legend('Est', 'True');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, ins_data(9,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, bacc_b_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Z acc bias [deg]')
legend('Est', 'True');
grid on;

% ATTITUDE 
figure(4)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, rad2deg*ins_data(10,:),'Color', 'blue', 'Linewidth', 2);
plot(time, rad2deg*att_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
ylabel('Roll angle [deg]')
legend('Est', 'True');
title('Attitude');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, rad2deg*ins_data(11,:), 'Color', 'blue', 'Linewidth', 2);
plot(time, rad2deg*att_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Pitch angle [deg]')
legend('Est', 'True');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, rad2deg*ins_data(12,:), 'Color', 'blue', 'Linewidth', 2);
plot(time, rad2deg*att_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
xlabel('Time [s]');
ylabel('yaw angle [deg]')
legend('Est', 'True');
grid on;

% GYRO BIAS
figure(5)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, ins_data(13,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, bars_b_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Roll bias [deg]')
legend('Est', 'True');
title('Gyro bias');
grid on;

subplot(3, 1, 2)
hold on;
plot(time, ins_data(14,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, bars_b_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Pitch bias [deg]')
legend('Est', 'True');
grid on;

subplot(3, 1, 3)
hold on;
plot(time, ins_data(15,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, bars_b_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Yaw bias [deg]')
legend('Est', 'True');
grid on;

% Position map
figure(6)
figure(gcf)
subplot(1, 1, 1)
plot(ins_data(2,:), ins_data(1,:), 'Color', 'black', 'Linewidth', 1.5);
xlabel('Y position [m]');
ylabel('X position [m]');
title('Position Plot');
grid on;
