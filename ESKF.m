clear all;
clc;
close all;

SKID = LoadVCUData('C:\Users\tonja\OneDrive\Skrivebord\CSV\FSG - Skid pad second stint',2,22883)

deg2rad = pi/180;
rad2deg = 180/pi;

simtime = 60;
f_samp = 200;
h = 1/f_samp;
N = simtime/h;

Z3 = zeros(3,3);
I3 = eye(3);

g = 9.81;

acc_bias = 2*[-0.4 -0.5 0.3]';
ars_bias = 10*[-0.03 0.02 -0.02]';

% INIT VALUES
x = [ones(1,3) zeros(1,6)];     % True state
xnoise = x;                     
delta_x = zeros(15,1);        % Error state
x_ins = zeros(15,1);          % Nominal state
xdot_ins = zeros(15,1);       

std_pos = 2;
std_att = 5 * deg2rad;

std_acc = 0.05;
std_acc_bias = 0.005;
std_ars = 0.04;
std_ars_bias = 0.005;


% KALMAN FILTER INIT
ESKF(0,0,0,10,false);

% ALLOCATE DATA STORAGE
time_data = zeros(1,N);
sim_data = zeros(


% --------------------------------------------------------------------- %
%                               FUNCTIONS                               %
% --------------------------------------------------------------------- %

function delta_x = ESKF(delta_y, R_nb, Tt, f_low, isInitialized)

    deg2rad = pi/180;
    
    Z3 = zeros(3,3);
    I3 = eye(3);
    
    Tacc = 200;
    Tars = 200;
    
    persistent P_hat Q R
    
    if (~isInitialized)
        std_pos = 2;
        R_pos = std_pos^2 * I3; 
        std_att = 10 * deg2rad;
        R_att = std_att^2 * I3;
        R = blkdiag(R_pos, R_att);
        
        std_acc = 0.05;
        Q_acc = ars_acc^2 * I3;
        std_acc_bias = 0.001;
        Q_acc_bias = std_acc_bias^2 * I3;
        
        std_ars = 0.1;
        Q_ars = ars_acc^2 * I3;
        std_ars_bias = 0.005;
        Q_ars_bias = std_acc_bias^2 * I3;
        
        % Q = blkdiag(Q_acc, Q_acc_bias, Q_ars, Q_ars_bias);
        Q = diag([1e-5 * ones(1, 3) 1e-2 * ones(1, 3) 1e-4 * ones(1, 3) 1e-5 * ones(1, 3)]);
        
        P_hat = diag([1e-6 * ones(1, 3) 1e-4 * ones(1, 3) 1e-4 * ones(1, 3) 1e-3 * ones(1, 3) 1e-3 * ones(1, 3)]);  % Initial error covariance
    
        delta_x = 0;
        
    else % not init
        
        % Error dynamics
        
        A = [ Z3  I3       Z3  Z3       Z3 ;
              Z3  Z3 -R_nb*I3  Z3       Z3 ;
              Z3  Z3 -I3/Tacc  Z3       Z3 ;
              Z3  Z3       Z3  Z3      -Tt ;
              Z3  Z3       Z3  Z3 -T3/Tars ];
          
        E = [  Z3  Z3  Z3  Z3 ;
              -I3  Z3  Z3  Z3 ;
               Z3  I3  Z3  Z3 ;  
               Z3  Z3 -I3  Z3 ;
               Z3  Z3  Z3  I3 ];
           
        C = [ I3  Z3  Z3  Z3  Z3 ;
              Z3  Z3  Z3  I3  Z3 ];
          
        % Discret time model
        [Ad, Ed] = c2d(A, E, 1/f_low);
        
        % Kalman gain
        K = P_hat * C' / (C * P_hat * C' + R);
        
        % Corrector
        delta_x = K * delta_y;
        P_hat = (eye(15) - K*C) * P_hat * (eye(15) - K*C)' + K*R*K';
        P_hat = (P_hat + P_hat') / 2;
        
        % Covariance predictor (k+1)
        P_hat = Ad * P_hat * Ad' + Ed * Q * Ed;
    
    end
          
           

end