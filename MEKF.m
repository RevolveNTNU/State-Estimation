clear all;
clc;
close all;

%% INITIALIZATION 
deg2rad = pi/180;   
rad2deg = 180/pi;
g       = [0 0 9.81]';

simtime = 100;
h       = 1/50;         %50 Hz update
N       = simtime/h;    %number of iterations

Z3 = zeros(3,3);
I3 = eye(3,3);

% MASS SPRING DAMPER
k = 1; d = 1; m = 1;
gm = 10; r = 1; l = 1;

%Plant initial values
phi = 0;
theta = 0;
psi = 0;
%init rotation matrix
[~, R_nb, Tt] = eulerang(phi, theta, psi);

x = [ones(1, 3) zeros(1, 3) phi theta psi zeros(1, 3)]';
x_est = zeros(12, 1);
x_ins = x;
xdot_ins = zeros(12, 1);

%init quaternion
q_ins = euler2q(phi, theta, psi);
omega_b_imu = zeros(3, 1);
f_b_imu = -g;

pos_meas = [0 0 0]';
newMeas = false;
measDelay = 0.2;

Tars = 1e4;

std_pos = .5;
std_att = 1 * deg2rad;
std_acc = 0.05;
std_ars = 0.04;
ars_bias = [-.030 0.02 -.04]';

R = diag([1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-1, 1e-1, 1e-1]);
Q = diag([1e-2 * ones(1, 3) 1e-4 * ones(1, 3) 1e-5 * ones(1, 3)]);

%initialization of Covariance matrix
P_hat = diag([1e-8 * ones(1, 3) 1e-8 * ones(1, 3) 1e-4 * ones(1, 3) 1e-4 * ones(1, 3)]);

time_data = zeros(1, N);
sim_data = zeros(12, N);
bias_data = zeros(3, N);
ins_data = zeros(12, N);
angvel_est = zeros(3, N);

HzCount10 = -measDelay/h;

%% MAIN LOOP
for i = 1:N
    t = (i-1) * h;                     
    ars_bias = ars_bias + 0.0005 * wgn(3, 1, 1);
    
    % Store data in simdata
    time_data(i) = t;
    sim_data(:,i) = x;
    ins_data(:,i) = x_ins;
    bias_data(:,i) = rad2deg * ars_bias; %acc_bias not used
    
    [phi, theta, psi] = q2euler(q_ins/norm(q_ins));
    ins_data(7:9, i) = [phi theta psi]';
    
    %INS measurements
    p_ins = x_ins(1:3);
    v_ins = x_ins(4:6);
    eps_ins = x_ins(7:9);
    q_ins = qbuild(eps_ins);
    q_ins    = q_ins/norm(q_ins);
    b_ins_ars = x_ins(10:12);
    R_ins = Rquat(q_ins);
    [~, ~, psi_ins] = q2euler(q_ins/norm(q_ins));
    
    %% plant

    %Plant rotation matrix
    phi_model     = x(7);
    theta_model   = x(8);
    psi_model     = x(9);
    [J, R_nb_model, Tt_model] = eulerang(phi_model, theta_model, psi_model);
    
    %model
    A_model = [Z3       I3      Z3          Z3      %p
            -k/m*I3     -d/m*I3     Z3      Z3      %v
            Z3          Z3          Z3      Tt_model*I3   %attitude
            Z3          Z3         -gm/l*I3 -d*I3];  %angular velocities

    B_model =   [Z3             Z3
                R_nb_model*I3/m       Z3 
                Z3              Z3
                Z3              Tt_model*I3/(m*r)];
            
    C_model  = [I3 Z3 Z3 Z3];
    
    %update 
    u = [5*[1 0.8 1.2]' * sin(.2*t) - g; 2 * [1 -0.8 1.1]' * sin(.5 * t)];
    xdot = A_model * x + B_model * u;
    
    a_b_nb  = xdot(4:6);
    acc_noise = std_acc*randn(3, 1);
    
    ars_noise = std_ars*randn(3, 1);
    
    omega_b_imu  = Tt_model' * x(10:12);
    f_b_imu = a_b_nb + skew(omega_b_imu) * x(4:6) + acc_noise + R_nb_model' * -g;
    omega_b_imu  = Tt_model' * x(10:12) + ars_bias + ars_noise;
    %% Error model
    A =    [Z3 I3 Z3 Z3 
            Z3 Z3 -R_ins * skew(f_b_imu) Z3
            Z3 Z3 -0.5 * skew(omega_b_imu - b_ins_ars) -0.5*I3
            Z3 Z3 Z3 -I3/Tars];

    E =    [ Z3 Z3 Z3
            -R_ins Z3 Z3 
             Z3 -I3 Z3
             Z3 Z3 I3];

    C = [I3 Z3 Z3 Z3
        Z3 I3 Z3 Z3
        Z3 Z3 I3 Z3];

    % Discrete-time model
    [Ad, Ed] = c2d(A, E, 0.1);   
    
    %% Get position measurements (10 Hz)
    
    %fast simulator used for GNSS delayed measurements
    A_fs = [I3 Z3; Z3 Z3];
    B_fs = [Z3 I3]';
    C_fs = [I3 Z3 Z3 Z3
            Z3 I3 Z3 Z3];
    x_fs = C_fs * x;
    
    % 10 Hz measurement update
    HzCount10 = HzCount10 + 1;
    if HzCount10 >= 1/h/10
        newMeas = true;
        HzCount10 = 0;
    end
    
    if (newMeas == true)          
        %delayed measurements (0.2 sec)
        tau = measDelay;        
        x_delayed = sim_data(:, i-tau/h);
        x_fs = C_fs * x_delayed;

        %integration to predict current position. 
        for j = 1:0.2/h
            x_fs_dot = A_fs * x_fs + B_fs * (f_n_imu + g);
            x_fs = x_fs + h * x_fs_dot;
        end
        
        pos_meas = x_fs(1:3) + std_pos * randn(1) * ones(3, 1);
        vel_meas = x_fs(4:6) + 0.001 * rand(1) * ones(3, 1);
       
        %% q-method:
     
        %reference vectors
        c_b = [cos(psi_ins)     sin(psi_ins)        0]';
        c_n = [1                0                   0]';
        
        q_met = qmethod(1,1,f_n_imu/norm(f_n_imu),c_n/norm(c_n) ,f_b_imu/norm(f_b_imu),c_b/norm(c_b));
        if norm(q_met+q_ins) < norm(q_met-q_ins) 
            q_met = -q_met;
        end
        
        %Find error from measurement
        dq = qmult(qinv(q_ins),q_met);
     
        dy = [pos_meas - p_ins; 
            vel_meas - v_ins; 
            dq(2:4)];

        %% KF update
        % KF gain
        K = P_hat * C' / (C * P_hat * C' + R);
        % new x_hat
        x_est = K * dy;
        % -----------
        p_est = x_est(1:3);
        v_est = x_est(4:6);
        eps_est = x_est(7:9);
        q_est = qbuild(eps_est);
        b_ins_ars_hat = x_est(10:12);
        % -----------
        
        % Covariance update 
        P_hat = (eye(12)-K*C) * P_hat * (eye(12) - K*C)' + K*R*K';
        P_hat = (P_hat + P_hat')/2;
        
        %Move error
        p_ins = p_ins + p_est;
        v_ins = v_ins + v_est;
        q_ins = qmult(q_ins, q_est);
        q_ins = q_ins/norm(q_ins);
        b_ins_ars = b_ins_ars + b_ins_ars_hat;

        x_ins = [p_ins; v_ins; q_ins(2:4); b_ins_ars];

        % Reset:
        x_est = zeros(12, 1);   

        % Covariance predictor (k+1)
        P_hat = Ad * P_hat * Ad' + Ed * Q * Ed';
        
        newMeas = false;
    end
    
    f_n_imu = R_ins * f_b_imu;
    
    %%
    % Strapdown INS equations
    omega_b_imu = omega_b_imu;
    xdot_ins(1:3) = x_ins(4:6);
    xdot_ins(4:6) = f_n_imu + g;
    qins_dot = 0.5 * qmult(q_ins, [0; omega_b_imu - b_ins_ars]);
    xdot_ins(7:9) = qins_dot(2:4);
    xdot_ins(10:12) = 0;
    
    x_fs_dot = A_fs * x_fs + B_fs * (f_n_imu + g);
    
    % Euler integration (k+1)
    x = x + h * xdot;
    x_ins = x_ins + h * xdot_ins;
    x_fs = x_fs + h * x_fs_dot;
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
q_ins           = ang_velocity(2, :);
r               = ang_velocity(3, :);

b_ars           = bias_data(1:3, :);

%ins measurements
p_ins           = ins_data(1:3, :);
xpos_ins        = p_ins(1, :);
ypos_ins        = p_ins(2, :);
zpos_ins        = p_ins(3, :);
v_ins           = ins_data(4:6, :);
xvel_ins        = v_ins(1, :);
yvel_ins        = v_ins(2, :);
zvel_ins        = v_ins(3, :);
attitude_ins    = rad2deg * ins_data(7:9, :);
phi_ins         = attitude_ins(1, :);
theta_ins       = attitude_ins(2, :);
psi_ins         = attitude_ins(3, :);
b_ars_ins       = rad2deg * ins_data(10:12, :);

figure(1)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(t, xpos);
plot(t, xpos_ins);
ylabel('X position [m]')
legend('True', 'Est');
title('Position');

subplot(3, 1, 2)
hold on;
plot(t, ypos);
plot(t, ypos_ins);
ylabel('Y position [m]')
legend('True', 'Est');

subplot(3, 1, 3)
hold on;
plot(t, zpos);
plot(t, zpos_ins);
xlabel('Time [s]');
ylabel('Z position [m]')
legend('True', 'Est');

% ATTITUDE 
figure(2)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(t, phi);
plot(t, phi_ins);
ylabel('Roll angle [deg]')
legend('True', 'Est');
title('Attitude');

subplot(3, 1, 2)
hold on;
plot(t, theta);
plot(t, theta_ins);
ylabel('Pitch angle [deg]')
legend('True', 'Est');

subplot(3, 1, 3)
hold on;
plot(t, psi);
plot(t, psi_ins);
xlabel('Time [s]');
ylabel('yaw angle [deg]')
legend('True', 'Est');

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

% GYRO BIAS
figure(4)
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

% POSITION MAP
figure(6)
figure(gcf)
subplot(1, 1, 1)
plot(xpos, ypos, 'Color', 'black', 'Linewidth', 1.5);
%legend('bias', 'est bias');
%legend('True', 'Est');
xlabel('X position [m]');
ylabel('Y position [m]');
title('Position Plot');
%%
function S = skew(x)
S = [   0 -x(3)  x(2) ;
             x(3)     0 -x(1) ;
            -x(2)  x(1)     0 ];
end

function q = qmethod(a1, a2, r1, r2, b1, b2)
Z = a1*cross(r1,b1) + a2*cross(r2,b2);
B = a1*r1*b1' + a2*r2*b2';
S = B + B';
sigma = trace(B);

% Compute the K matrix
K = [ -sigma     Z'
       Z        -S+sigma*eye(3) ]; 

% Find the eigenvector for the smallest eigenvalue of K
[V,E] = eig(K);
[m,i] = min([E(1,1) E(2,2) E(3,3) E(4,4)]);
q = V(:,i);
q = q/norm(q);

end

function q = qmult(q1, q2)
if any(~isreal(q1(:)))
    error('First input elements are not real.');
end
if any(size(q1)~=[4,1])
    error('First input dimension is not 4-by-1.');
end
if any(~isreal(q2(:)))
    error('Second input elements are not real.');
end
if any(size(q2)~=[4,1])
    error('Second input dimension is not 4-by-1.');
end

eta1 = q1(1);
eps1 = q1(2:4);
eta2 = q2(1);
eps2 = q2(2:4);

eta = eta1*eta2 - eps1'*eps2;
eps = eta1*eps2 + eta2*eps1 + cross(eps1,eps2);

q = [eta; eps];

end

function qinv_out = qinv(q)
q(2:4) = -q(2:4);
qinv_out = q;
end

function q = qbuild(eps)
eta = sqrt(1 - eps'*eps);
q   = [eta; eps];
end




