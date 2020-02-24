clear all;
clc;
close all;

deg2rad = pi/180;   
rad2deg = 180/pi;

Z3 = zeros(3,3);
I3 = eye(3);

simtime = 300;
f_samp  = 100;          %imu frequency
h       = 1/f_samp;     %sampling time
N       = simtime/h;    %number of iterations

%init values
delta_x = zeros(15, 1);
x_ins = zeros(15, 1);
y = zeros(6,1);

%initialization of kalman filter
ErrorStateKalman(0, 0, 0, 1);

% data storage
time_data = zeros(1, N);
ins_data = zeros(15, N);

% from sim
[p_n_nb, v_n_nb, att_n_nb, f_b_imu, omega_b_imu,time, g_b_n, v_abs] = CircleSim(simtime,f_samp,0);

% init
x_ins(1:3) = [5;6;7]; % for testing av ESKF
x_ins(4:6) = [0;0;0]; % for testing av ESKF
x_ins(10:12) = [deg2rad*13; -deg2rad*15; -deg2rad*42]; % for testing av ESKF
% x_ins(4:6) = [v_abs;0;0]; % korrekt initialisering. For testing av INS

% matrices
C_ins = [I3 Z3 Z3 Z3 Z3
         Z3 Z3 Z3 I3 Z3];
     
     
count = 10;

for k = 1:N
   t = k * h;
   
   time_data(k) = t;
   ins_data(:,k) = x_ins;

   phi_ins     = x_ins(10);
   theta_ins   = x_ins(11);
   psi_ins     = x_ins(12);
    
   [J_ins, R_nb_ins, Tt_ins] = eulerang(phi_ins, theta_ins, psi_ins);
   
   count = count + 1;

   if (count >= 10)
        count = 0;
        
        % perfect measurement
        y(1:3) = p_n_nb(1:3,k);
        y(4:6) = att_n_nb(1:3,k);
        
        y_ins = C_ins * x_ins;
        delta_y = y - y_ins;
        
        delta_x = ErrorStateKalman(delta_y, R_nb_ins, Tt_ins, 0);
        x_ins = x_ins + delta_x;
        
   end
   

    
   a_n_ins = R_nb_ins*f_b_imu(:,k) + g_b_n;
   omega_b_ins = omega_b_imu(:,k) - x_ins(13:15);
 
   x_ins(1:3) = x_ins(1:3) + (h * x_ins(4:6)) + (0.5 * h * h * a_n_ins);
   x_ins(4:6) = x_ins(4:6) +  (h * a_n_ins);
   x_ins(7:9) = x_ins(7:9);
   x_ins(10:12) = x_ins(10:12) + (h * Tt_ins*omega_b_ins);
   x_ins(13:15) = x_ins(13:15); 
   
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

subplot(3, 1, 2)
hold on;
plot(time, ins_data(2,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, p_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Y position [m]')
legend('Est', 'True');

subplot(3, 1, 3)
hold on;
plot(time, ins_data(3,:) ,'Color', 'blue', 'Linewidth', 2);
plot(time, p_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
ylabel('Z position [m]')
legend('Est', 'True');

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


% ATTITUDE 
figure(3)
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


% Position map
figure(4)
figure(gcf)
subplot(1, 1, 1)
plot(ins_data(2,:), ins_data(1,:), 'Color', 'black', 'Linewidth', 1.5);
xlabel('Y position [m]');
ylabel('X position [m]');
title('Position Plot');
grid on;

