clear all;
clc;
close all;
% 
% rosbag_moholt_autox_2;
% rosbag_moholt_autox_good;

deg2rad = pi/180;   
rad2deg = 180/pi;

simtime = 240;
f_samp  = 100;          %imu frequency
h       = 1/f_samp;     %sampling time
N       = simtime/h;    %number of iterations

Z3 = zeros(3,3);
I3 = eye(3);

% MASS SPRING DAMPER
k = 1;
d = 1;
m = 1;

g = 10;
r = 1;
l = 1;

acc_bias = 2*[-0.4 -0.5 0.3]';
ars_bias = 10*[-.030 0.02 -.02]';

%init values
x = [ones(1, 3) zeros(1, 9)]';
xnoise = x;
delta_x = [zeros(15, 1)];
x_ins = [zeros(15, 1)];
xdot_ins = [zeros(15, 1)];

std_pos = 2;
std_att = 5 * deg2rad;

std_acc = 0.05;
std_acc_bias = 0.005;

std_ars = 0.04;
std_ars_bias = 0.005;

%initialization of kalman filter
IndirectKalman(0, 0, 0, 1);

time_data = zeros(1, N);
sim_data = zeros(12, N);
sim_data_noise = zeros(12, N);
bias_data = zeros(6, N);
ins_data = zeros(15, N);
angvel_est = zeros(3, N);


% Simulator

% choose only one case
% 1 - Standing still at position [100,100,0]' with attitide
%     deg2rad*[0,0,45]'
% 2 - Driving on a straight line starting at position [10,10,0]' with 
%     constant velocity of 10 m/s and orientation deg2rad*[0,0,45]
% 3 - Vehicle model provided in TK8109 based on mass-spring-samper system
%     with sinusoidal acceleration and angular rate input
% 4 - Circular motion 

SIM_CASE = 4;

if (SIM_CASE == 1)
    x(1:3) = [100 100 0]'; % position
    x(7:9) = [0 0 deg2rad*45]'; % attitude

elseif (SIM_CASE == 2)
    % Create position measurements
    pos = zeros(3,N);
    pos(1:3, 1) = [10 10 0]';

    for i = 2:N
        pos(1:3, i) = pos(1:3, i-1) + 0.01*[5*sqrt(2) 5*sqrt(2) 0]';
    end
    
    % Modify init values
    x(1:3) = [10 10 0]'; % position
    x(4:6) = [5*sqrt(2) 5*sqrt(2) 0]';
    x(7:9) = deg2rad*[0 0 45]';    
elseif (SIM_CASE == 3)
    % keep initi values
elseif (SIM_CASE == 4)
    R = 50;
    abs_vel = (2 * R * pi) / 30;

    pos = zeros(3,N);
    pos(1:3, 1) = [0,-50,0]';

    att = zeros(3,N);
    att(1:3, 1) = [0,0,0]';
    
    alpha = 0;
    for i = 2:N
        alpha = alpha + (pi/1500); 
        pos(1:3, i) = [R*sin(alpha) -R*cos(alpha) 0]';
        att(1:3, i) = [0 0 alpha]';
    end
    
    x(1:3) = [0 -50 0]';
    x(4:6) = [abs_vel 0 0]';
    x(7:9) = [0 0 0]';
    
    x_ins(1:3) = [0 -50 0]';
    x_ins(4:6) = [abs_vel 0 0]';
    x_ins(7:9) = [0 0 0]';
    
    
    
end



count = 10;
for i = 2:N
    t = (i-1) * h;   
    
    ars_bias = ars_bias + 0.00005 * wgn(3, 1, 1);
    acc_bias = acc_bias + 0.001 * wgn(3, 1, 1);
    % Store data in simdata
    time_data(i) = t;
    sim_data(:,i) = x;
    sim_data(4:6,i) = sim_data(4:6,i); 
    sim_data_noise(:,i) = xnoise;
    ins_data(:,i) = x_ins;
    %bias_data(:,i) = [acc_bias; rad2deg*ars_bias];
    
    %euler angles
    phi     = x(7);
    theta   = x(8);
    psi     = x(9);
    
    [J, R_nb, Tt] = eulerang(phi, theta, psi);
    
    
    
    
%     % Strapdown INS
%     A_ins = [Z3 I3 Z3 Z3 Z3         %p
%              Z3 Z3 -R_nb*I3 Z3 Z3   %v
%              Z3 Z3 Z3 Z3 Z3         %acc_bias
%              Z3 Z3 Z3 Z3 -Tt*I3     %theta
%              Z3 Z3 Z3 Z3 Z3];       %ars_bias
% 
%     B_ins = [Z3 Z3
%              R_nb*I3 Z3
%              Z3 Z3
%              Z3 Tt*I3
%              Z3 Z3];

     C_ins = [I3 Z3 Z3 Z3 Z3
              Z3 Z3 Z3 I3 Z3];
% 
%     %plant 
% 
%     if (SIM_CASE == 1 || SIM_CASE == 2) 
%     A_model = [ Z3  I3  Z3  Z3
%                 Z3  Z3  Z3  Z3
%                 Z3  Z3  Z3  I3
%                 Z3  Z3  Z3  Z3] ;
% 
% 
%     B_model =   [Z3             Z3
%                 R_nb*I3         Z3 
%                 Z3              Z3
%                 Z3              Tt*I3];
% 
%     C_model  = [I3 Z3 Z3 Z3
%                 Z3 Z3 I3 Z3];
% 
% 
%     elseif (SIM_CASE == 3 || SIM_CASE == 4 )        
%         % Mass-spring-damper
% 
%         A_model = [Z3       I3      Z3          Z3      %p
%                 -k/m*I3     -d/m*I3     Z3      Z3      %v
%                 Z3          Z3          Z3      Tt*I3   %attitude
%                 Z3          Z3         -g/l*I3 -d*I3];  %angular velocities
% 
%         B_model =   [Z3             Z3
%                     R_nb*I3/m       Z3 
%                     Z3              Z3
%                     Z3              Tt*I3/(m*r)];

        C_model  = [I3 Z3 Z3 Z3
                    Z3 Z3 I3 Z3];
%     end

            
    % 10 Hz update
    count = count + 1;
%     if count >= 10
%         count = 0;
%         
%         % Measurements
%         xnoise = x + [std_pos * randn(1) * ones(1, 3) 0 0 0 std_att * randn(1) * ones(1, 3) 0 0 0]';
%         %y = C_model * x + [std_pos * randn(1) * ones(1, 3) std_att * randn(1) * ones(1, 3)]';             
%         y_ins = C_ins*x_ins;
%         
%         if (SIM_CASE == 1)
%             y(1:3,1) = [100, 100, 0];
%             y(4:6,1) = [0, 0, deg2rad*45];
%         elseif (SIM_CASE == 2)
%             y(1:3,1) = pos(:,i);
%             y(4:6,1) = [0 , 0, 45]' * deg2rad;
%         elseif (SIM_CASE == 3)
%             y = C_model * x;
%         elseif (SIM_CASE == 4)
%             y(1:3,1) = pos(:,i);
%             y(4:6,1) = att(:,i);
%         end
%         
%         y = y + [std_pos * randn(1) * ones(1, 3), std_att * randn(1) * ones(1, 3)]';
%         
% %         [px,py,pz] = llh2flat(SKID.vcu_GNSS_longitude(i),SKID.vcu_GNSS_latitude(i),SKID.vcu_GNSS_altitude(i),SKID.vcu_GNSS_longitude(2),SKID.vcu_GNSS_latitude(2),SKID.vcu_GNSS_altitude(2)); 
%         
% %         y(1) = 0.2 * -px;
% %         y(2) = 0.2 * -py;
% %         y(3) = 0.2 * -pz;
% %         
% %         y(4) = deg2rad * SKID.vcu_INS_roll(i);
% %         y(5) = deg2rad * SKID.vcu_INS_pitch(i);
% %         y(6) = deg2rad * SKID.vcu_INS_yaw(i);
%          
%         
%         delta_y = y - y_ins;
%         
%         delta_x = IndirectKalman(delta_y, R_nb, Tt, 0);
%         
%         x_ins = x_ins + delta_x;
%     end
    
    % true input
    if (SIM_CASE == 1 || SIM_CASE == 2)
         u = [ [0 0 0]' ; [0 0 0]' ]; 
    elseif (SIM_CASE == 3 ) 
        u = [5*[1 0.8 1.2]' * sin(.2*t); 2 * [1 -0.8 1.1]' * sin(.5 * t)];
    elseif (SIM_CASE == 4)
        u = [0*[-(abs_vel^2)/R * sin(x(9)), (abs_vel^2)/R * cos(x(9)), 0]' ; [0,0,pi/15]']; 
    end
        
    % u = [SKID.vcu_INS_ax(i) SKID.vcu_INS_ay(i) SKID.vcu_INS_az(i) SKID.vcu_INS_roll_rate(i) SKID.vcu_INS_pitch_rate(i) SKID.vcu_INS_yaw_rate(i)]'; 
    
    %xdot = A_model * x + B_model * u;
    
    a_b_nb  = [0 0 0]';
    acc_noise = std_acc*randn(3, 1);
    w_b_nb  = [0 0 pi/15]';
    ars_noise = std_ars*randn(3, 1);
       
    
    % f_imu_b = [SKID.vcu_INS_ax(i) SKID.vcu_INS_ay(i) SKID.vcu_INS_az(i)]';
    % w_imu = [SKID.vcu_INS_roll_rate(i) SKID.vcu_INS_pitch_rate(i) SKID.vcu_INS_yaw_rate(i)]';
    
    % Strapdown INS equations
     %xdot_ins = A_ins * x_ins + B_ins * [f_imu_b;(w_b_nb + ars_bias + ars_noise)];
     
    
%     xdot_ins(1:3) = x_ins(4:6);
%     xdot_ins(4:6) = R_nb*(f_imu_b - acc_bias) + [0 0 g]';
%     xdot_ins(7:9) = zeros(1,3)';
%     xdot_ins(10:12) = Tt*(w_b_nb + ars_noise);
%     xdot_ins(13:15) = zeros(1,3)';
    
    f_b_imu = a_b_nb + ScrewSym(u(4:6)) * (R_nb') * x(4:6); % + acc_bias + acc_noise;
    w_b_imu = w_b_nb ; %+ ars_bias + ars_noise;
    
    
    phi_ins     = x_ins(7);
    theta_ins   = x_ins(8);
    psi_ins     = x_ins(9);
    
    [J_ins, R_nb_ins, Tt_ins] = eulerang(phi_ins, theta_ins, psi_ins);
    

    a_n_ins = R_nb_ins*(f_b_imu - acc_bias);
    
    x_ins(1:3) = x_ins(1:3) + (h * x_ins(4:6)) + (0.5 * h * h * a_n_ins);
    x_ins(4:6) = x_ins(4:6) +  (h * a_n_ins);
    x_ins(7:9) = x_ins(7:9);
    x_ins(10:12) = x_ins(10:12) + h* Tt_ins *(w_b_imu - ars_bias);
    x_ins(13:15) = x_ins(13:15); 
    
    
         
    x(1:3) = x(1:3) + (h * x(4:6));
    x(4:6) = x(4:6) + (h * R_nb * u(1:3));
    x(7:9) = x(7:9) + (h * Tt * u(4:6));
    % x(10:12) = ((1 - h)* x(10:12)) - (h * x(7:9)) + (h * Tt * u(4:6)); 
    
    % Euler integration (k+1)
    %x = x + h * xdot;
    %x_ins = x_ins + h * xdot_ins;
    
    
end

%% PLOTS
t               = time_data;

p               = sim_data(1:3, :);
xpos            = p(1, :);
ypos            = p(2, :);
zpos            = p(3, :);

v               = sim_data(4:6, :);
xvel            = v(1, :);
yvel            = v(2, :);
zvel            = v(3, :);

attitude        = rad2deg * sim_data(7:9, :);
phi             = attitude(1, :);
theta           = attitude(2, :);
psi             = attitude(3, :);

ang_velocity    = sim_data(10:12, :);
p               = ang_velocity(1, :);
q               = ang_velocity(2, :);
r               = ang_velocity(3, :);

b_acc           = bias_data(1:3, :);
b_ars           = bias_data(4:6, :);

% noisy measurements 
p_n             = sim_data_noise(1:3, :);
xpos_n          = p_n(1, :);
ypos_n          = p_n(2, :);
zpos_n          = p_n(3, :);

attitude_n      = rad2deg * sim_data_noise(7:9, :);
phi_n           = attitude_n(1, :);
theta_n         = attitude_n(2, :);
psi_n           = attitude_n(3, :);

%ins measurements
p_ins           = ins_data(1:3, :);
xpos_ins        = p_ins(1, :);
ypos_ins        = p_ins(2, :);
zpos_ins        = p_ins(3, :);

v_ins           = ins_data(4:6, :);
xvel_ins        = v_ins(1, :);
yvel_ins        = v_ins(2, :);
zvel_ins        = v_ins(3, :);

b_acc_ins       = ins_data(7:9, :);

attitude_ins    = rad2deg * ins_data(10:12, :);
phi_ins         = attitude_ins(1, :);
theta_ins       = attitude_ins(2, :);
psi_ins         = attitude_ins(3, :);

b_ars_ins       = rad2deg*ins_data(13:15, :);

p_ins           = angvel_est(1, :);
q_ins           = angvel_est(2, :);
r_ins           = angvel_est(3, :);

% POSITION
figure(1)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(t, xpos, 'Color', 'black', 'Linewidth', 1.5);
plot(t, xpos_n, 'Color', [255/255,165/255,0]);
plot(t, xpos_ins, '--', 'Color', 'blue', 'Linewidth', 2);
ylabel('X position [m]')
legend('True', 'Meas', 'Est');
title('Position');

subplot(3, 1, 2)
hold on;
plot(t, ypos, 'Color', 'black', 'Linewidth', 1.5);
plot(t, ypos_n, 'Color', [255/255,165/255,0]);
plot(t, ypos_ins, '--', 'Color', 'blue', 'Linewidth', 2);
ylabel('Y position [m]')
legend('True', 'Meas', 'Est');

subplot(3, 1, 3)
hold on;
plot(t, zpos, 'Color', 'black', 'Linewidth', 1.5);
plot(t, zpos_n, 'Color', [255/255,165/255,0]);
plot(t, zpos_ins, '--', 'Color', 'blue', 'Linewidth', 2);
xlabel('Time [s]');
ylabel('Z position [m]')
legend('True', 'Meas', 'Est');


% ATTITUDE 
figure(2)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(t, phi, 'Color', 'black', 'Linewidth', 1.5);
plot(t, phi_n, 'Color', [255/255,165/255,0]);
plot(t, phi_ins, '--', 'Color', 'blue', 'Linewidth', 2);
ylabel('Roll angle [deg]')
legend('True', 'Meas', 'Est');
title('Attitude');

subplot(3, 1, 2)
hold on;
plot(t, theta, 'Color', 'black', 'Linewidth', 1.5);
plot(t, theta_n, 'Color', [255/255,165/255,0]);
plot(t, theta_ins, '--', 'Color', 'blue', 'Linewidth', 2);
ylabel('Pitch angle [deg]')
legend('True', 'Meas', 'Est');

subplot(3, 1, 3)
hold on;
plot(t, psi, 'Color', 'black', 'Linewidth', 1.5);
plot(t, psi_n, 'Color', [255/255,165/255,0]);
plot(t, psi_ins, '--', 'Color', 'blue', 'Linewidth', 2);
xlabel('Time [s]');
ylabel('yaw angle [deg]')
legend('True', 'Meas', 'Est');

% VELOCITIES
figure(3)
figure(gcf)
subplot(3, 1, 1)
hold on;
plot(t, xvel);
plot(t, xvel_ins);
xlabel('Time [s]');
ylabel('X velocity [m/s]')
legend('True', 'Est');
title('Velcoity');

subplot(3, 1, 2)
hold on;
plot(t, yvel);
plot(t, yvel_ins);
legend('True', 'Est');

subplot(3, 1, 3)
hold on;
plot(t, zvel);
plot(t, zvel_ins);
legend('True', 'Est');

% ACCELERATION BIAS
figure(4)
figure(gcf)
subplot(3, 1, 1)
plot(t, b_acc(1, :), t, b_acc_ins(1, :));
legend('True', 'Est');
ylabel('Xacc bias [m/s^2]')
title('Accelerometer bias');

subplot(3, 1, 2)
plot(t, b_acc(2, :), t, b_acc_ins(2, :));
legend('True', 'Est');
ylabel('Yacc bias [m/s^2]')

subplot(3, 1, 3)
plot(t, b_acc(3, :), t, b_acc_ins(3, :));
legend('True', 'Est');
xlabel('Time [s]');
ylabel('Zacc bias [m/s^2]')

% GYRO BIAS
figure(5)
figure(gcf)
subplot(3, 1, 1)
plot(t, b_ars(1, :), t, b_ars_ins(1, :));
%legend('bias', 'est bias');
legend('True', 'Est');
ylabel('Roll bias [deg]')
title('Gyro bias');

subplot(3, 1, 2)
plot(t, b_ars(2, :), t, b_ars_ins(2, :));
%legend('bias', 'est bias');
legend('True', 'Est');
ylabel('Pitch bias [deg]')

subplot(3, 1, 3)
plot(t, b_ars(3, :), t, b_ars_ins(3, :));
legend('True', 'Est');
xlabel('Time [s]');
ylabel('Yaw bias [deg]')

if (SIM_CASE == 2)
    % POSITION MAP
    figure(6)
    figure(gcf)
    subplot(1, 1, 1)
    plot(pos(2,:), pos(1,:), 'Color', 'black', 'Linewidth', 1.5);
    %legend('bias', 'est bias');
    %legend('True', 'Est');
    xlabel('Y position [m]');
    ylabel('X position [m]');
    title('Position Plot');

% 
%     % ATTITUDE MAP
%     figure(7)
%     figure(gcf)
%     subplot(1, 1, 1)
%     plot(t, att(3,:), 'Color', 'black', 'Linewidth', 1.5);
%     %legend('bias', 'est bias');
%     %legend('True', 'Est');
%     xlabel('time [s]');
%     ylabel('yaw [m]');
%     title('Attitude Plot');

end

if (SIM_CASE == 4)
    % POSITION MAP
    figure(6)
    figure(gcf)
    subplot(1, 1, 1)
    plot(-ypos_ins, xpos_ins, 'Color', 'black', 'Linewidth', 1.5);
    xlabel('Y position [m]');
    ylabel('X position [m]');
    title('Position Plot');
end
    
%% functions

function S = ScrewSym(x)
S = [   0 -x(3)  x(2) ;
             x(3)     0 -x(1) ;
            -x(2)  x(1)     0 ];
end


