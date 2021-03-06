function [delta_x, E] = ErrorStateKalman_sola(f_b_ins, race_started, r_b_1, r_b_2, r_b_3, E_prev,delta_y, R_nb_hat, f_low, init, f_b_imu, omega_b_imu, g_n_nb, x_n_ins)
    
    deg2rad = pi/180;   
    
    Z3 = zeros(3,3);
    I3 = eye(3); 
    
    Tacc = 3600;
    Tars = 3600;
    
    h = 1/f_low; 
    
    p_n_ins = x_n_ins(1:3);
    v_n_ins = x_n_ins(4:6);
    bacc_b_ins = x_n_ins(7:9);
    q_n_ins = x_n_ins(10:13);
    q_n_ins = q_n_ins/norm(q_n_ins);
    bars_b_ins = x_n_ins(14:16);
        
    persistent P_hat Q R  
    if init
        std_pos = 2;
        R_pos = std_pos^2*I3;
        std_att = 1 * deg2rad;
        R_att = std_att^2*I3;
        std_vel = 1;
        R_vel = std_vel^2*I3;
        std_acc = 1;
        R_acc = std_acc^2*I3;
        R = blkdiag(R_pos , R_pos, std_vel^2, 2*R_pos, R_acc);
%         R = blkdiag(R_pos , R_pos);


        std_acc = 0.01 * sqrt(10);
        %std_acc = 0.14 * 0.001 * 9.80665;
        Q_acc = std_acc * std_acc * I3; 
        std_acc_bias = 1;
        %std_acc_bias = 0.04 * 0.001 * 9.80665;
        Q_acc_bias = std_acc_bias * std_acc_bias * I3;

        std_ars = 0.1;
        %std_ars = 0.0035 * deg2rad;
        Q_ars = std_ars * std_ars * I3;
        std_ars_bias = 0.01 * sqrt(10);
        %std_ars_bias = 10 * deg2rad;
        Q_ars_bias = std_ars_bias * std_ars_bias * I3;

        Q = blkdiag( Q_acc, Q_acc_bias, Q_ars, Q_ars_bias ) * h * h;
        Q = 100 * Q;
        
%         std_acc = 0.01 * sqrt(10);
%         std_acc = 0.14 * 0.001 * 9.80665;
%         Q_acc = std_acc * std_acc * h *I3; 
%         std_acc_bias = 1;
%         std_acc_bias = 0.04 * 0.001 * 9.80665;
%         Q_acc_bias = 2 * std_acc_bias * std_acc_bias * h * I3;
% 
%         std_ars = 0.1;
%         std_ars = 0.0035 * deg2rad;
%         Q_ars = std_ars * std_ars * h * I3;
%         std_ars_bias = 0.01 * sqrt(10);
%         std_ars_bias = 10 * deg2rad;
%         Q_ars_bias = 2 * std_ars_bias * std_ars_bias * h * I3;
% 
%         Q = blkdiag( Q_acc, Q_acc_bias, Q_ars, Q_ars_bias ); % * h * h;
%        
        
%         Q = diag([1e-5 * ones(1, 3) 1e-2 * ones(1, 3) 1e-4 * ones(1, 3) 1e-5 * ones(1, 3)]);
      
%        P_hat = diag([1e-6 * ones(1, 3) 1e-4 * ones(1, 3) 1e-4 * ones(1, 3) 1e-3 * ones(1, 3) 1e-3 * ones(1, 3)]);  % Initial error covariance

%         P_hat = diag([ 1e-1 1e-1 1e-1 ...
%                        1e-2 1e-2 1e-2 ...
%                        5e-2 5e-2 5e-2 ...
%                        1e-10 1e-10 1e-10 ...
%                        1e-6 1e-6 1e-6]); 
                   
       P_hat = diag([1e-1 * ones(1, 3) 1e-2 * ones(1, 3) 5e-2 * ones(1, 3) 1e-10 * ones(1, 3) 1e-6 * ones(1, 3) 1e-6 * ones(1, 3)]);  % Initial error covariance

       delta_x = 0;
       
    else
          A_dp = [Z3  I3  Z3  Z3  Z3  Z3];
          A_dv = [Z3  Z3  -R_nb_hat  -R_nb_hat*Smtrx(f_b_imu - bacc_b_ins)  Z3  I3];
          A_dbacc = [Z3  Z3  -(1/Tacc)*I3  Z3  Z3  Z3];
          A_dtheta = [Z3  Z3  Z3  -Smtrx(omega_b_imu - bars_b_ins)  -I3  Z3];
          A_dbars = [Z3  Z3  Z3  Z3  -(1/Tars)*I3  Z3];
          A_dg = [Z3  Z3  Z3  Z3  Z3  Z3];
          
          if (~race_started)
              % Do not estimate 
              A_dv = [Z3  Z3  -R_nb_hat  -R_nb_hat*Smtrx(f_b_imu - bacc_b_ins)  Z3  Z3];
          end
          
          A =  [A_dp ; A_dv ; A_dbacc ; A_dtheta ; A_dbars ; A_dg];
          
          
         % Discrete-time model
         Ad = eye(18) + h * A;
%          if (~race_started)
%              % Do not estimate gravity error
%              Ad = blkdiag(I3,I3,I3,I3,I3,Z3) + h * A;
%          end
          
%           A = [  Z3    I3            Z3                                       Z3           Z3  Z3   % dp
%                  Z3    Z3     -R_nb_hat    -R_nb_hat*Smtrx(f_b_imu - bacc_b_ins)           Z3  I3   % dv
%                  Z3    Z3  -(1/Tacc)*I3                                       Z3           Z3  Z3   % dbacc
%                  Z3    Z3            Z3         -Smtrx(omega_b_imu - bars_b_ins)          -I3  Z3   % dtheta
%                  Z3    Z3            Z3                                       Z3 -(1/Tars)*I3  Z3   % dbars
%                  Z3    Z3            Z3                                       Z3           Z3  Z3]; % dg
%              
            
 
          E = [        Z3      Z3    Z3   Z3
                -R_nb_hat      Z3    Z3   Z3    % w_acc
                       Z3      I3    Z3   Z3    % w_acc_bias
                       Z3      Z3   -I3   Z3    % w_ars
                       Z3      Z3    Z3   I3  % w_ars_bias
                       Z3      Z3    Z3   Z3 ]; 

          
%           H = [ I3 Z3 Z3 Z3 Z3 Z3
%                 I3 Z3 Z3 I3 Z3 Z3];

          H_gnss1 = [I3  Z3  Z3  -R_nb_hat*Smtrx(r_b_1)  Z3  Z3];
          
          H_gnss2 = [I3  Z3  Z3  -R_nb_hat*Smtrx(r_b_2)  Z3  Z3];
          
          H_vec = [Z3  Z3  Z3  -Smtrx(R_nb_hat*(r_b_2-r_b_1))  Z3  Z3];
          
          H_acc = [Z3  Z3  Z3  -Smtrx(R_nb_hat' * g_n_nb)  Z3  Z3];
%           H_acc = [Z3  Z3  I3  -Smtrx(R_nb_hat' * g_n_nb)  Z3  Z3];   
%           H_acc = [Z3  Z3  Z3  Z3  Z3  Z3]; 

          H_gss_pos = [ 0, 0, 0];
          H_gss_vel = [(R_nb_hat(1,2)*(R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) + R_nb_hat(1,1)*(R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), (R_nb_hat(2,2)*(R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) + R_nb_hat(1,3)*(R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), (R_nb_hat(3,2)*(R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) + R_nb_hat(2,3)*(R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2)];
          H_gss_bacc = [0, 0, 0]; 
          H_gss_att = [((R_nb_hat(1,3)*v_n_ins(1) + R_nb_hat(2,3)*v_n_ins(2) + R_nb_hat(3,3)*v_n_ins(3))*(R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)))/((R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -((R_nb_hat(1,3)*v_n_ins(1) + R_nb_hat(2,3)*v_n_ins(2) + R_nb_hat(3,3)*v_n_ins(3))*(R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -(2*(R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3))*(R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) - 2*(R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3))*(R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/(2*((R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2))];
          H_gss_bars = [(r_b_3(3)*(R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)))/((R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -(r_b_3(3)*(R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -(r_b_3(1)*(R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) - r_b_3(2)*(R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_hat(1,2)*v_n_ins(1) + R_nb_hat(2,2)*v_n_ins(2) + R_nb_hat(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_hat(1,1)*v_n_ins(1) + R_nb_hat(1,3)*v_n_ins(2) + R_nb_hat(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2)]; 
          H_gss_g = [0, 0, 0];
          
          H_gss = [H_gss_pos  H_gss_vel  H_gss_bacc  H_gss_att  H_gss_bars  H_gss_g];
          
          
        std_pos = 2;
        R_pos1 = std_pos^2*I3;
        R_pos2 = std_pos^2*I3;
        R_vec = 20*std_pos^2*I3;
        std_att = 1 * deg2rad;
        R_att = std_att^2*I3;
        std_vel = 1;
        R_vel = std_vel^2;
        std_acc = 1;
        R_acc = std_acc^2*I3;
        R = blkdiag(R_pos1 , R_pos2, std_vel^2, R_vec, R_acc);
          
           % Outlier Rejection
          [~,H_gnss1] = ChiSquareTest(H_gnss1, P_hat, R_pos1, delta_y(1:3));
          [~,H_gnss2] = ChiSquareTest(H_gnss2, P_hat, R_pos2, delta_y(4:6));
          [~,H_gss] = ChiSquareTest(H_gss, P_hat, R_vel, delta_y(7));
          [~,H_vec] = ChiSquareTest(H_vec, P_hat, R_vec, delta_y(8:10));
          [~,H_acc] = ChiSquareTest(H_acc, P_hat, R_acc, delta_y(11:13)); 
          
          
          if (~race_started)
            H = [H_gnss1 ; H_gnss2 ; H_vec ; H_acc];
            R = blkdiag(R_pos1 , R_pos2, R_vec, R_acc);
            delta_y(7) = []; % delete gss 
%             H_gss = zeros(1,18);
%             H_acc = [Z3  Z3  I3  -Smtrx(R_nb_hat' * g_n_nb)  Z3  -R_nb_hat'];
%             H_acc = [Z3  Z3  I3  -Smtrx(R_nb_hat' * g_n_nb)  Z3  Z3];    
          else
            H = [H_gnss1 ; H_gnss2 ; H_gss ; H_vec];
            R = blkdiag(R_pos1 , R_pos2, R_vel, R_vec);
            delta_y(11:13) = []; % delete acc 
          end 
          
    
%           
%           H = [        I3         Z3          Z3        -R_nb_hat*Smtrx(r_b_1)          Z3          Z3   % gnss_1
%                        I3         Z3          Z3        -R_nb_hat*Smtrx(r_b_2)          Z3          Z3   % gnss_2
%                 H_gss_pos  H_gss_vel  H_gss_bacc                     H_gss_att  H_gss_bars     H_gss_g % ground speed
% %                        Z3         Z3          Z3                      I3          Z3       Z3]; % delta_theta
%                        Z3         Z3          Z3  -R_nb_hat*Smtrx(r_b_2-r_b_1)          Z3          Z3   % baseline
%                        Z3         Z3          I3     -Smtrx(R_nb_hat' * g_n_nb)          Z3  -R_nb_hat']; % acc
%                    
%           if (race_started == true)
%              H = H(1:end-3, :);
%              if (size(R,1) == 13)
%                 R = R(1:end-3, 1:end-3);
%              end
%           end
              
               
%          M = rref(obsv(A,H))  
         r = rank(obsv(A,H));
%          if (~race_started)
%             disp(r); 
%          end
          

        
         % KF gain
         K = P_hat * H' / (H * P_hat * H' + R);

         % corrector 
         delta_x = K * delta_y;
         P_hat = (eye(18)-K*H) * P_hat * (eye(18) - K*H)' + K*R*K';
        
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
        P_hat = (P_hat + P_hat')/2;

             
    end
end