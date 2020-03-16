function [delta_x, E] = ErrorStateKalman2(theta_der,Q_deltaq, E_prev,delta_y, R_nb, f_low, init, f_b_imu, omega_b_imu, g_n_nb, bacc_b_ins, bars_b_ins)
    deg2rad = pi/180;   
    
    Z3 = zeros(3,3);
    I3 = eye(3);
    
    Tacc = 200;
    Tars = 200;
    
    h = 1/f_low;        

    persistent P_hat Q R  
    if init
        std_pos = 2;
        R_pos = std_pos^2*I3;
        std_att = 10 * deg2rad;
        R_att = std_att^2*I3;
        R = blkdiag(R_pos, R_att);

        std_acc = 0.01 * sqrt(10);
        Q_acc = std_acc * std_acc * I3; 
        std_acc_bias = 1;
        Q_acc_bias = std_acc_bias * std_acc_bias * I3;

        std_ars = 0.1;
        Q_ars = std_ars * std_ars *I3;
        std_ars_bias = 0.01 * sqrt(10);
        Q_ars_bias = std_ars_bias * std_ars_bias * I3;

        Q = blkdiag( Q_acc, Q_acc_bias, Q_ars, Q_ars_bias ) * h * h;
%         Q = diag([1e-5 * ones(1, 3) 1e-2 * ones(1, 3) 1e-4 * ones(1, 3) 1e-5 * ones(1, 3)]);
      
%        P_hat = diag([1e-6 * ones(1, 3) 1e-4 * ones(1, 3) 1e-4 * ones(1, 3) 1e-3 * ones(1, 3) 1e-3 * ones(1, 3)]);  % Initial error covariance

%         P_hat = diag([ 1e-1 1e-1 1e-1 ...
%                        1e-2 1e-2 1e-2 ...
%                        5e-2 5e-2 5e-2 ...
%                        1e-10 1e-10 1e-10 ...
%                        1e-6 1e-6 1e-6]); 
                   
       P_hat = diag([1e-1 * ones(1, 3) 1e-2 * ones(1, 3) 5e-2 * ones(1, 3) 1e-10 * ones(1, 3) 1e-6 * ones(1, 3) 1e-10 * ones(1, 3)]);  % Initial error covariance

       delta_x = 0;
%        delta_q = 0;
       
    else
          omega_b_ins = omega_b_imu - bars_b_ins;
          phi = h * omega_b_ins(1);
          theta = h * omega_b_ins(2);
          psi = h * omega_b_ins(3);
          R_omega = Rzyx(phi, theta, psi);

            
          A = [  Z3    I3        Z3                                   Z3       Z3    Z3    % dp
                 Z3    Z3     -R_nb    -R_nb*Smtrx(f_b_imu - bacc_b_ins)       Z3    I3    % dv
                 Z3    Z3        Z3                                   Z3       Z3    Z3    % dbacc
                 Z3    Z3        Z3      -Smtrx(omega_b_imu - bars_b_ins)     -I3    Z3    % dtheta ???????
                 Z3    Z3        Z3                                   Z3       Z3    Z3    % dbars
                 Z3    Z3        Z3                                   Z3       Z3    Z3] ; % dg
            
 
          E = [      Z3      Z3    Z3   Z3
                  -R_nb      Z3    Z3   Z3    % w_acc
                     Z3      I3    Z3   Z3    % w_acc_bias
                     Z3      Z3   -I3   Z3    % w_ars
                     Z3      Z3    Z3   I3    % w_ars_bias
                     Z3      Z3    Z3   Z3 ]; 
                 
          
             
          X = blkdiag(I3,I3,I3,Q_deltaq,I3,I3);
          Q_34 = [ 0 1 0 0
                   0 0 1 0
                   0 0 0 1];
          Z_34 = zeros(3,4);
               
%           H = [ I3 Z3 Z3 Z_34 Z3 Z3
%                 Z3 Z3 Z3 Q_34 Z3 Z3];
            
          
          H = [ I3 Z3 Z3 zeros(3,4) Z3 Z3
                Z3 Z3 Z3  theta_der Z3 Z3];
            
          C = H * X;
          
          
%           C = [ I3 Z3 Z3 Z3 Z3 Z3];
          
          % Discrete-time model
          Ad = eye(18) + h * A;
%           Ed = h * E;
        
        % KF gain
        K = P_hat * C' / (C * P_hat * C' + R);

        % corrector 
        delta_x = K * delta_y;
        P_hat = (eye(18)-K*C) * P_hat * (eye(18) - K*C)' + K*R*K';
        P_hat = (P_hat + P_hat')/2;
        
        % ESKF reset
        delta_theta = delta_x(10:12);
        G = [ I3 Z3 Z3                            Z3 Z3 Z3 
              Z3 I3 Z3                            Z3 Z3 Z3 
              Z3 Z3 I3                            Z3 Z3 Z3 
              Z3 Z3 Z3 (I3 - Smtrx(0.5*delta_theta)) Z3 Z3 
              Z3 Z3 Z3                            Z3 I3 Z3 
              Z3 Z3 Z3                            Z3 Z3 I3];
          
        P_hat = G * P_hat * G';
        
        
        % Covariance predictor (k+1)

        Qd = 0.5 * (Ad * E_prev * Q * E_prev' * Ad' + E * Q * E') * h;
        P_hat = Ad * P_hat * Ad' + Qd;
        
%         P_hat = Ad * P_hat * Ad' + Ed * Q * Ed';

     
    end
end