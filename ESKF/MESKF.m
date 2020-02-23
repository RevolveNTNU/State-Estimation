clear all;
clc;
close all;

deg2rad = pi/180;   
rad2deg = 180/pi;

Z3 = zeros(3,3);
I3 = eye(3);

simtime = 240;
f_samp  = 100;          %imu frequency
f_low   = 10;           %aiding frequency
h       = 1/f_samp;     %sampling time
N       = simtime/h;    %number of iterations

% init values
delta_x = zeros(15, 1);
x_ins = zeros(15, 1);
y = zeros(6,1);

% init quaternion
phi_0 = deg2rad*45;
theta_0 = deg2rad*45;
psi_0 = deg2rad*45;
q_ins = euler2q(phi_0,theta_0,psi_0); % q = [eta eps1 eps2 eps3]

% init rotation matrix
%[~, R_nb, Tt] = eulerang(phi_0, theta_0, psi_0);
R_nb = Rquat(q_ins);
Tt = Tquat(q_ins);

std_pos = 2;
std_att = 5 * deg2rad;

%initialization of kalman filter
f_b_imu_0 = [0 0 0];
w_b_imu_0 = [0 0 0];
ErrorStateKalman2(0, 0, 0, f_low, 1,f_b_imu_0 , w_b_imu_0, q_ins);

% data storage
time_data = zeros(1, N);
ins_data = zeros(15, N);

% from sim
[p_n_nb, v_n_nb, att_n_nb, f_b_imu, w_b_imu,time] = CircleSim(simtime,f_samp,0);
% [p_n_nb, v_n_nb, att_n_nb, f_b_imu, w_b_imu,time] = StandStillSim(simtime,f_samp,0);


% matrices
C_ins = [I3 Z3 Z3 Z3 Z3
         Z3 Z3 Z3 I3 Z3];
     
     
count = 10;

for i = 2:N
   t = i * h;
   
   time_data(i) = t;
   ins_data(:,i) = x_ins;

   [phi_ins, theta_ins, psi_ins] = q2euler(q_ins/norm(q_ins)); 
%    phi_ins = ssa(phi_ins, 'rad');
%    theta_ins = ssa(theta_ins, 'rad');
%    psi_ins = ssa(psi_ins, 'rad');
   ins_data(7:9, i) = [phi_ins theta_ins psi_ins]'; %for plotting
   
%    phi_ins     = x_ins(7);
%    theta_ins   = x_ins(8);
%    psi_ins     = x_ins(9);
    
    p_n_ins = x_ins(1:3);
    v_n_ins = x_ins(4:6);
    acc_bias_ins = x_ins(7:9);
    eps_ins = x_ins(10:12);
    ars_bias_ins = x_ins(13:15);
    
    % eps_ins = x_ins(10:12);    
    q_ins = unitquat(eps_ins);
    q_ins = q_ins/norm(q_ins); % normalization to avoid numerical erro

   [J_ins, R_nb_ins, Tt_ins] = eulerang(phi_ins, theta_ins, psi_ins);

   R_nb_ins = Rquat(q_ins);
   Tt_ins = Tquat(q_ins);
   
   count = count + 1;

   if (count >= 10)
        count = 0;
        
        % perfect measurement
        y(1:3) = p_n_nb(1:3,i); % + (std_pos * randn(1) * ones(1, 3))';
        y(4:6) = att_n_nb(1:3,i); % + (std_pos * randn(1) * ones(1, 3))';
        
        y_ins = C_ins * x_ins;
        delta_y = y - y_ins;
        
        [delta_x, delta_q] = ErrorStateKalman2(delta_y, R_nb_ins, Tt_ins, f_low, 0, f_b_imu, w_b_imu, q_ins);
        % x_ins = x_ins + delta_x;
        
        
        p_n_ins = p_n_ins + delta_x(1:3);
        v_n_ins = v_n_ins + delta_x(4:6);
        acc_bias_ins = acc_bias_ins + delta_x(7:9);
        q_ins = quatprod(q_ins, delta_q);
        q_ins = q_ins / norm(q_ins);
        ars_bias_ins = ars_bias_ins + delta_x(13:15);
        
        
   end
   
   a_n_ins = R_nb_ins*f_b_imu;
   
   % from Sola kinematics
  

    p_n_ins = p_n_ins + (h * v_n_ins) + (0.5 * h * h * a_n_ins);
    v_n_ins = v_n_ins + (h * a_n_ins);
    acc_bias_ins = acc_bias_ins;
    q_ins = q_ins + 0.5 * h * quatprod(q_ins, [0 ; w_b_imu]); % egt (w_b_imu - ars_bias)
    ars_bias_ins = ars_bias_ins; 
    
    x_ins = [p_n_ins; v_n_ins; acc_bias_ins; q_ins(2:4); ars_bias_ins];
    
%     x_ins(1:3) = x_ins(1:3) + (h * x_ins(4:6)) + (0.5 * h * h * a_n_ins);
%     x_ins(4:6) = x_ins(4:6) +  (h * a_n_ins);
%     x_ins(7:9) = x_ins(7:9);
%     x_ins(10:12) = x_ins(10:12) + (h * w_b_imu);
%     x_ins(13:15) = x_ins(13:15); 
    
end


   % PLOTS
   
% POSITION
figure(1)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time_data, ins_data(1,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time_data, p_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
ylabel('X position [m]')
legend('Est', 'True');
title('Position');

subplot(3, 1, 2)
hold on;
plot(time_data, ins_data(2,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time_data, p_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Y position [m]')
legend('Est', 'True');

subplot(3, 1, 3)
hold on;
plot(time_data, ins_data(3,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time_data, p_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Z position [m]')
legend('Est', 'True');

% VELOCITIES
figure(2)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time_data, ins_data(4,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time_data, v_n_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('X velocity [m]')
legend('Est', 'True');
title('Velocity');

subplot(3, 1, 2)
hold on;
plot(time_data, ins_data(5,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time_data, v_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Y velocity [m]')
legend('Est', 'True');

subplot(3, 1, 3)
hold on;
plot(time_data, ins_data(6,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time_data, v_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Z velocity [m]')
legend('Est', 'True');

% ATTITUDE 
figure(3)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time_data, rad2deg*ins_data(10,:),'Color', 'blue', 'Linewidth', 2);
plot(time_data, rad2deg*att_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
ylabel('Roll angle [deg]')
legend('Est', 'True');
title('Attitude');

subplot(3, 1, 2)
hold on;
plot(time_data, rad2deg*ins_data(11,:), 'Color', 'blue', 'Linewidth', 2);
plot(time_data, rad2deg*att_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Pitch angle [deg]')
legend('Est', 'True');

subplot(3, 1, 3)
hold on;
plot(time_data, rad2deg*ins_data(12,:), 'Color', 'blue', 'Linewidth', 2);
plot(time_data, rad2deg*att_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
xlabel('Time [s]');
ylabel('yaw angle [deg]')
legend('Est', 'True');

% Position map
figure(4)
figure(gcf)
subplot(1, 1, 1)
plot(-ins_data(2,:), ins_data(1,:), 'Color', 'black', 'Linewidth', 1.5);
xlabel('Y position [m]');
ylabel('X position [m]');
title('Position Plot');

