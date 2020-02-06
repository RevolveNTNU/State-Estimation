% clear all;
% clc;
% close all;
% 
% rosbag_moholt_autox_2;
% rosbag_moholt_autox_good;

deg2rad = pi/180;   
rad2deg = 180/pi;

simtime = 90;
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

%i_gps = 29000;
l0 = NOVA.vcu_GNSS_longitude(1);
mu0 = NOVA.vcu_GNSS_latitude(1);
h_ref = NOVA.vcu_GNSS_altitude(1);

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

count = 10;
for i = 1:N
    t = (i-1) * h;   
    
    ars_bias = ars_bias + 0.00005 * wgn(3, 1, 1);
    acc_bias = acc_bias + 0.001 * wgn(3, 1, 1);
    % Store data in simdata
    time_data(i) = t;
    sim_data(:,i) = x;
    sim_data_noise(:,i) = xnoise;
    ins_data(:,i) = x_ins;
    bias_data(:,i) = [acc_bias; rad2deg*ars_bias];
    
    %euler angles
    phi     = x(7);
    theta   = x(8);
    psi     = x(9);
    [J, R_nb, Tt] = eulerang(phi, theta, psi);

    % Strapdown INS
    A_ins = [Z3 I3 Z3 Z3 Z3         %p
             Z3 Z3 -R_nb*I3 Z3 Z3   %v
             Z3 Z3 Z3 Z3 Z3         %acc_bias
             Z3 Z3 Z3 Z3 -Tt*I3     %theta
             Z3 Z3 Z3 Z3 Z3];       %ars_bias
         
    B_ins = [Z3 Z3
             R_nb*I3 Z3
             Z3 Z3
             Z3 Tt*I3
             Z3 Z3];
         
     C_ins = [I3 Z3 Z3 Z3 Z3
              Z3 Z3 Z3 I3 Z3];

    %plant
    A_model = [Z3       I3      Z3          Z3      %p
            -k/m*I3     -d/m*I3     Z3      Z3      %v
            Z3          Z3          Z3      Tt*I3   %attitude
            Z3          Z3         -g/l*I3 -d*I3];  %angular velocities

    B_model =   [Z3             Z3
                R_nb*I3/m       Z3 
                Z3              Z3
                Z3              Tt*I3/(m*r)];

    C_model  = [I3 Z3 Z3 Z3
                Z3 Z3 I3 Z3];
            
    % 10 Hz update
    count = count + 1;
    if count >= 10
        count = 0;
        %i_gps = i_gps + 1;
        
        % Measurements
        xnoise = x + [std_pos * randn(1) * ones(1, 3) 0 0 0 std_att * randn(1) * ones(1, 3) 0 0 0]';
        y = C_model * x + [std_pos * randn(1) * ones(1, 3) std_att * randn(1) * ones(1, 3)]';             
        y_ins = C_ins*x_ins;
        
        % GPS
        [px,py,pz] = llh2flat(NOVA.vcu_GNSS_longitude(i),NOVA.vcu_GNSS_latitude(i),NOVA.vcu_GNSS_altitude(i),l0,mu0,h_ref); 
        
        y(1) = 0,1*px + std_pos * randn(1);
        y(2) = 0.1*py + std_pos * randn(1);
        y(3) = 0.1*pz + std_pos * randn(1);
        
        delta_y = y - y_ins;
        disp(delta_y);
        
        delta_x = IndirectKalman(delta_y, R_nb, Tt, 0);
        disp(delta_x);
        
        x_ins = x_ins + delta_x;
    end
    
    %plant
    % u = [5*[1 0.8 1.2]' * sin(.2*t); 2 * [1 -0.8 1.1]' * sin(.5 * t)];
    u = [NOVA.vcu_INS_ax(i) NOVA.vcu_INS_ay(i) NOVA.vcu_INS_az(i) NOVA.vcu_INS_roll_rate(i) NOVA.vcu_INS_pitch_rate(i) NOVA.vcu_INS_yaw_rate(i)]'; 
    xdot = A_model * x + B_model * u;
    
    a_b_nb  = xdot(4:6);
    acc_noise = std_acc*randn(3, 1);
    w_b_nb  = xdot(7:9);
    ars_noise = std_ars*randn(3, 1);
       
    f_imu_b = a_b_nb + ScrewSym(x(10:12)) * x(7:9) + acc_bias + acc_noise;
    
    % Strapdown INS equations
    xdot_ins = A_ins * x_ins + B_ins * [f_imu_b; (w_b_nb + ars_bias + ars_noise)];
         
    % Euler integration (k+1)
    x = x + h * xdot;
    x_ins = x_ins + h * xdot_ins;
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
yvel_ins        = v_ins(2, :)   ;
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

%% functions

function S = ScrewSym(x)
S = [   0 -x(3)  x(2) ;
             x(3)     0 -x(1) ;
            -x(2)  x(1)     0 ];
end

function delta_x = IndirectKalman(delta_y, R_nb, Tt, init)
    deg2rad = pi/180;   
    
    Z3 = zeros(3,3);
    I3 = eye(3);
    
    Tacc = 200;
    Tars = 200;

    persistent P_hat Q R
    if init
        std_pos = 2;
        R_pos = std_pos^2*I3;
        std_att = 10 * deg2rad;
        R_att = std_att^2*I3;
        R = blkdiag(R_pos, R_att);

        std_acc = 0.05;
        Q_acc = std_acc^2*I3; 
        std_acc_bias = 0.001;
        Q_acc_bias = std_acc_bias^2*I3;

        std_ars = 0.1;
        Q_ars = std_ars^2*I3;
        std_ars_bias = 0.005;
        Q_ars_bias = std_ars_bias^2*I3;

        Q = blkdiag( Q_acc, Q_acc_bias, Q_ars, Q_ars_bias );
        Q = diag([1e-5 * ones(1, 3) 1e-2 * ones(1, 3) 1e-4 * ones(1, 3) 1e-5 * ones(1, 3)]);
      
        P_hat = diag([1e-6 * ones(1, 3) 1e-4 * ones(1, 3) 1e-4 * ones(1, 3) 1e-3 * ones(1, 3) 1e-3 * ones(1, 3)]);  % Initial error covariance

        delta_x = 0;
    else
            % Error model
        A =    [Z3 I3  Z3 Z3 Z3 
                Z3 Z3 -R_nb*I3 Z3 Z3
                Z3 Z3 -I3/Tacc Z3 Z3
                Z3 Z3 Z3 Z3 -Tt
                Z3 Z3 Z3 Z3 -I3/Tars];

        E =    [ Z3 Z3 Z3 Z3
                -I3 Z3 Z3 Z3 
                 Z3 I3 Z3 Z3
                 Z3 Z3 -I3 Z3
                 Z3 Z3 Z3 I3];

        C = [I3 Z3 Z3 Z3 Z3
            Z3 Z3 Z3 I3 Z3];

        % Discrete-time model
        [Ad, Ed] = c2d(A, E, 0.1);

        % KF gain
        K = P_hat * C' / (C * P_hat * C' + R);

        % corrector 
        delta_x = K * delta_y;
        P_hat = (eye(15)-K*C) * P_hat * (eye(15) - K*C)' + K*R*K';
        P_hat = (P_hat + P_hat')/2;

        % Covariance predictor (k+1)
        P_hat = Ad * P_hat * Ad' + Ed * Q * Ed';
    end
end
