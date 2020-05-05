clear all;
clc;
close all;

deg2rad = pi/180;   
rad2deg = 180/pi;

Z3 = zeros(3,3);
I3 = eye(3);

simtime = 250;
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
[p_n_nb, v_n_nb, q_nb, bacc_b_nb, bars_b_nb, f_b_imu, omega_b_imu,time, g_n_nb, v_abs, acc_std, ars_std] = SkidPadSim(simtime,f_samp,0);
% [p_n_nb, v_n_nb, att_n_nb, f_b_imu, w_b_imu,time] = StandStillSim(simtime,f_samp,0);


% dual gnss config
r_b_1 = [0.5 0.5 0.5]';
r_b_2 = [-0.5 -0.2 0]';
r_b_3 = [0.5 -0.4 0.1]'; 

%initialization of kalman filter
f_b_imu_0 = [0 0 0]';
omega_b_imu_0 = [0 0 0]';
bacc_b_ins = [0 0 0]';
bars_b_ins = [0 0 0]';
E_prev = zeros(18,12);
ErrorStateKalman_sola(0,0,r_b_1, r_b_2, r_b_3, E_prev,0, 0, f_low, 1, f_b_imu_0, omega_b_imu_0, g_n_nb, x_ins);
Iterative_ESKF(0, 0, 0, 0, r_b_1, r_b_2, r_b_3, E_prev,0, 0, f_low, 1, f_b_imu_0,omega_b_imu_0, g_n_nb, x_ins);
% init
x_ins(1:3) = [5;6;7]; % for testing av ESKF
% x_ins(1:3) = p_n_nb(:,1);
x_ins(4:6) = [0;0;0]; % for testing av ESKF
% x_ins(4:6) = v_n_nb(:,1);
x_ins(10:13) = [1 0 0 0]'; % for testing av ESKF
% x_ins(4:6) = [v_abs;0;0]; % korrekt initialisering. For testing av INS

% matrices
C_ins = [I3 Z3 Z3 Z3 Z3
         Z3 Z3 Z3 I3 Z3];
     
g_n_hat = g_n_nb;
     
count = 10;
race_started = false;

for k = 1:N
   t = k * h;
   time_data(k) = t;
   g_err_data(1:3,k) = g_n_hat;
   
   if ((p_n_nb(1,k) - p_n_nb(1,1)) > 0.001) && (race_started == false)
       race_started = true;
       disp(t);
   end
   
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
   a_n_ins = R_nb_ins*(f_b_imu(:,k) - bacc_b_ins) + g_n_hat;
   f_b_ins = R_nb_ins' * (a_n_ins - g_n_hat);
   omega_b_ins = omega_b_imu(:,k) - bars_b_ins;
   
   % compute quaternion from angular rates 
   q_omega_b_ins = qbuild(omega_b_ins, h);
   
   % update nominal states with imu input
   p_n_ins = p_n_ins + (h * v_n_ins) + (0.5 * h * h * a_n_ins);
   v_n_ins = v_n_ins +  (h * a_n_ins);
  %bacc_b_ins = bacc_b_ins;
   q_n_ins = quatprod(q_n_ins, q_omega_b_ins); 
   q_n_ins = q_n_ins/norm(q_n_ins);
  %bars_b_ins = bars_b_ins;
   
   x_ins = [p_n_ins ; v_n_ins ; bacc_b_ins ; q_n_ins ; bars_b_ins];
   
   
   count = count + 1;

   if (count >= 10)
        count = 0;
        
        % noisy measurements
        p_meas(1:3,k) = p_n_nb(1:3,k) +  0.001 * randn(3, 1);
        q_meas = q_nb(1:4,k) + (2 * deg2rad) * randn(4, 1);
        
        q_conj = quatconj(q_n_ins')';
        delta_q = quatprod(q_conj, q_meas);
        delta_theta = 2*delta_q(2:4);
        
%         p_gnss_1 = p_n_nb(1:3,k) + I3*R_nb_ins*r_b_1 - Smtrx(R_nb_ins*r_b_1)*delta_theta;
%         p_gnss_2 = p_n_nb(1:3,k) + I3*R_nb_ins*r_b_2 - Smtrx(R_nb_ins*r_b_2)*delta_theta;

% -----------------------------------------------------------------------
        % Dual gnss position meas
        R_nb_t = Rquat(q_nb(1:4,k)');
        y_gnss1 = p_n_nb(1:3,k) + R_nb_t * r_b_1;
        y_gnss2 = p_n_nb(1:3,k) + R_nb_t * r_b_2;

        y_gnss1_hat = p_n_ins + R_nb_ins * r_b_1;
        y_gnss2_hat = p_n_ins + R_nb_ins * r_b_2;
        
        H_gnss1 = [I3  Z3  Z3  -R_nb_ins*Smtrx(r_b_1)  Z3  Z3];
        H_gnss2 = [I3  Z3  Z3  -R_nb_ins*Smtrx(r_b_2)  Z3  Z3];

        std_pos = 2;
        R_pos = std_pos^2*I3;
        
        [delta_x, E_prev] = Iterative_ESKF(H_gnss1, R_pos, f_b_ins, race_started, r_b_1, r_b_2, r_b_3, E_prev,(y_gnss1 - y_gnss1_hat), R_nb_ins, f_low, 0, f_b_imu, omega_b_imu, g_n_nb, x_ins);
        [x_ins,g_n_hat] = InjectErrorState(delta_x, x_ins, g_n_hat);
        
        [delta_x, E_prev] = Iterative_ESKF(H_gnss2, R_pos, f_b_ins, race_started, r_b_1, r_b_2, r_b_3, E_prev,(y_gnss2 - y_gnss2_hat), R_nb_ins, f_low, 0, f_b_imu, omega_b_imu, g_n_nb, x_ins);
        [x_ins,g_n_hat] = InjectErrorState(delta_x, x_ins, g_n_hat);
        

        
% -----------------------------------------------------------------------
%         % Ground speed velocity meas
        if (race_started)
            H_gss_alloc = [1 0 0; 0 1 0; 0 0 0];
            y_gss = H_gss_alloc*( R_nb_ins'*v_n_nb(1:3,k) + Smtrx( omega_b_imu(:,k) - bars_b_nb)*r_b_3 );
            y_gss = norm( y_gss );

            y_gss_hat = R_nb_ins' * v_n_ins - r_b_3(3)*(bars_b_ins(2) - omega_b_imu(2,k)) + r_b_3(2)*(bars_b_ins(3) - omega_b_imu(3,k));
            y_gss_hat = norm(y_gss_hat);

            H_gss_pos = [ 0, 0, 0];
            H_gss_vel = [(R_nb_ins(1,2)*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) + R_nb_ins(1,1)*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), (R_nb_ins(2,2)*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) + R_nb_ins(1,3)*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), (R_nb_ins(3,2)*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) + R_nb_ins(2,3)*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2)];
            H_gss_bacc = [0, 0, 0]; 
            H_gss_att = [((R_nb_ins(1,3)*v_n_ins(1) + R_nb_ins(2,3)*v_n_ins(2) + R_nb_ins(3,3)*v_n_ins(3))*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -((R_nb_ins(1,3)*v_n_ins(1) + R_nb_ins(2,3)*v_n_ins(2) + R_nb_ins(3,3)*v_n_ins(3))*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -(2*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3))*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) - 2*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3))*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/(2*((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2))];
            H_gss_bars = [(r_b_3(3)*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -(r_b_3(3)*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -(r_b_3(1)*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) - r_b_3(2)*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2)]; 
            H_gss_g = [0, 0, 0];

            H_gss = [H_gss_pos  H_gss_vel  H_gss_bacc  H_gss_att  H_gss_bars  H_gss_g];

            std_vel = 1;
            R_vel = std_vel^2;

            [delta_x, E_prev] = Iterative_ESKF(H_gss, R_vel, f_b_ins, race_started, r_b_1, r_b_2, r_b_3, E_prev,(y_gss - y_gss_hat), R_nb_ins, f_low, 0, f_b_imu, omega_b_imu, g_n_nb, x_ins);
            [x_ins, g_n_hat] = InjectErrorState(delta_x, x_ins, g_n_hat);
        end
        
        
% -----------------------------------------------------------------------
        % Moving baseline
        y_vec = R_nb_t * (r_b_2 - r_b_1) + 1*randn(1)*[5,1,1]';
        y_vec_hat = R_nb_ins * (r_b_2 - r_b_1);
        
        H_vec = [Z3  Z3  Z3  -Smtrx(R_nb_ins*(r_b_2-r_b_1))  Z3  Z3];
        
        std_pos = 2;
        R_vec = 20*std_pos^2*I3;

        [delta_x, E_prev] = Iterative_ESKF(H_vec, R_vec, f_b_ins, race_started, r_b_1, r_b_2, r_b_3, E_prev,(y_vec - y_vec_hat), R_nb_ins, f_low, 0, f_b_imu, omega_b_imu, g_n_nb, x_ins);
        [x_ins, g_n_hat] = InjectErrorState(delta_x, x_ins, g_n_hat);
        
        
% -----------------------------------------------------------------------        
%         % Stand-still acceleration
        if (~race_started)
            y_acc = f_b_imu(:,k); %-(R_nb_t)' * g_n_nb + bacc_b_nb(:,k);
            y_acc_hat = -(R_nb_ins)' * g_n_hat; % + bacc_b_ins;

            H_acc = [Z3  Z3  Z3  -Smtrx(R_nb_ins' * g_n_nb)  Z3  Z3];        

            std_acc = 1;
            R_acc = std_acc^2*I3;

            [delta_x, E_prev] = Iterative_ESKF(H_acc, R_acc, f_b_ins, race_started, r_b_1, r_b_2, r_b_3, E_prev,(y_acc - y_acc_hat), R_nb_ins, f_low, 0, f_b_imu, omega_b_imu, g_n_nb, x_ins);
            [x_ins, g_n_hat] = InjectErrorState(delta_x, x_ins, g_n_hat);
        end
        
        
        
%         % compute difference 
% %         delta_y = [(p_gnss_1 - p_hat_1) ; (p_gnss_2 - p_hat_2); (p_gnss_3 - p_hat_3); delta_theta]; 
% %         delta_y = [(p_gnss_1 - p_hat_1) ; (p_gnss_2 - p_hat_2); delta_theta]; 
% %         delta_y = [(p_gnss_1 - p_hat_1); (p_gnss_2 - p_hat_2)];
%         delta_y = [(y_gnss1 - y_gnss1_hat) ; (y_gnss2 - y_gnss2_hat) ; (y_gss - y_gss_hat) ; (y_vec - y_vec_hat) ; (y_acc - y_acc_hat)] ;
%         
%         % GNSS1
%         [delta_x, E_prev] = ErrorStateKalman_sola(f_b_ins, race_started, r_b_1, r_b_2,r_b_3, E_prev, delta_y(1:3), R_nb_ins, f_low, 0, f_b_imu(:,k), omega_b_imu(:,k), g_n_hat, x_ins);
%         x_ins = InjectErrorState(delta_x, x_ins, g_n_hat);
%         
%         % GNSS2
%         [delta_x, E_prev] = ErrorStateKalman_sola(f_b_ins, race_started, r_b_1, r_b_2,r_b_3, E_prev, delta_y(4:6), R_nb_ins, f_low, 0, f_b_imu(:,k), omega_b_imu(:,k), g_n_hat, x_ins);
%         x_ins = InjectErrorState(delta_x, x_ins, g_n_hat);
%         
%         % GSS
%         if (race_started)
%             [delta_x, E_prev] = ErrorStateKalman_sola(f_b_ins, race_started, r_b_1, r_b_2,r_b_3, E_prev, delta_y(7), R_nb_ins, f_low, 0, f_b_imu(:,k), omega_b_imu(:,k), g_n_hat, x_ins);
%             x_ins = InjectErrorState(delta_x, x_ins, g_n_hat);
%         end
%         
%         % VEC
%         [delta_x, E_prev] = ErrorStateKalman_sola(f_b_ins, race_started, r_b_1, r_b_2,r_b_3, E_prev, delta_y(8:10), R_nb_ins, f_low, 0, f_b_imu(:,k), omega_b_imu(:,k), g_n_hat, x_ins);
%         x_ins = InjectErrorState(delta_x, x_ins, g_n_hat);
%         
%         % ACC
%         if (~race_started)
%             [delta_x, E_prev] = ErrorStateKalman_sola(f_b_ins, race_started, r_b_1, r_b_2,r_b_3, E_prev, delta_y(11:13), R_nb_ins, f_low, 0, f_b_imu(:,k), omega_b_imu(:,k), g_n_hat, x_ins);
%             x_ins = InjectErrorState(delta_x, x_ins, g_n_hat);
%         end
%         
%         
%         
%         
%                   H_gnss1 = [I3  Z3  Z3  -R_nb_ins*Smtrx(r_b_1)  Z3  Z3];
%           
%           H_gnss2 = [I3  Z3  Z3  -R_nb_ins*Smtrx(r_b_2)  Z3  Z3];
%           
%           H_vec = [Z3  Z3  Z3  -Smtrx(R_nb_ins*(r_b_2-r_b_1))  Z3  Z3];
%           
%           H_acc = [Z3  Z3  Z3  -Smtrx(R_nb_ins' * g_n_nb)  Z3  Z3];              
% %           H_acc = [Z3  Z3  Z3  Z3  Z3  Z3]; 
% 
%           H_gss_pos = [ 0, 0, 0];
%           H_gss_vel = [(R_nb_ins(1,2)*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) + R_nb_ins(1,1)*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), (R_nb_ins(2,2)*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) + R_nb_ins(1,3)*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), (R_nb_ins(3,2)*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) + R_nb_ins(2,3)*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2)];
%           H_gss_bacc = [0, 0, 0]; 
%           H_gss_att = [((R_nb_ins(1,3)*v_n_ins(1) + R_nb_ins(2,3)*v_n_ins(2) + R_nb_ins(3,3)*v_n_ins(3))*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -((R_nb_ins(1,3)*v_n_ins(1) + R_nb_ins(2,3)*v_n_ins(2) + R_nb_ins(3,3)*v_n_ins(3))*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -(2*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3))*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) - 2*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3))*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/(2*((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2))];
%           H_gss_bars = [(r_b_3(3)*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -(r_b_3(3)*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2), -(r_b_3(1)*(R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1)) - r_b_3(2)*(R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2)))/((R_nb_ins(1,2)*v_n_ins(1) + R_nb_ins(2,2)*v_n_ins(2) + R_nb_ins(3,2)*v_n_ins(3) - bars_b_ins(2)*r_b_3(1) + bars_b_ins(1)*r_b_3(3) - omega_b_imu(1)*r_b_3(3) + omega_b_imu(3)*r_b_3(1))^2 + (R_nb_ins(1,1)*v_n_ins(1) + R_nb_ins(1,3)*v_n_ins(2) + R_nb_ins(2,3)*v_n_ins(3) + bars_b_ins(2)*r_b_3(2) - bars_b_ins(2)*r_b_3(3) + omega_b_imu(2)*r_b_3(3) - omega_b_imu(3)*r_b_3(2))^2)^(1/2)]; 
%           H_gss_g = [0, 0, 0];
%           
%           H_gss = [H_gss_pos  H_gss_vel  H_gss_bacc  H_gss_att  H_gss_bars  H_gss_g];
%           
%           
%         std_pos = 2;
%         R_pos1 = std_pos^2*I3;
%         R_pos2 = std_pos^2*I3;
%         R_vec = 2*std_pos^2*I3;
%         std_att = 1 * deg2rad;
%         R_att = std_att^2*I3;
%         std_vel = 1;
%         R_vel = std_vel^2;
%         std_acc = 1;
%         R_acc = std_acc^2*I3;
%         R = blkdiag(R_pos1 , R_pos2, std_vel^2, R_vec, R_acc);
%         
%         
%         
%         
%         
%         
%         
%         
%         
%         
%         
        
        
        
%         delta_y = [(y_gnss1 - y_gnss1_hat) ; (y_gnss2 - y_gnss2_hat) ; (y_gss - y_gss_hat) ; (y_vec - y_vec_hat)];
        
%         if (race_started == false)
%             delta_y = [delta_y ; (f_imu - f_hat)];
%         end
        
%         f_b_ins = Smtrx(omega_b_imu(:,k) - bars_b_ins) * R_nb_ins' * v_n_ins;
        
        % compute error state with ESKF
%         [delta_x, E_prev] = ErrorStateKalman_sola(f_b_ins, race_started, r_b_1, r_b_2,r_b_3, E_prev, delta_y, R_nb_ins, f_low, 0, f_b_imu(:,k), omega_b_imu(:,k), g_n_hat, x_ins);
       
%         disp(delta_x);
        
%         % inject error state into nominal state
%         x_ins(1:9) = x_ins(1:9) + delta_x(1:9);
%         x_ins(14:16) = x_ins(14:16) + delta_x(13:15);
%         g_n_hat = g_n_hat + delta_x(16:18);
% %         disp(delta_x(16:18));
%         h_low = 1/10;
%         q_delta_omega = qbuild(delta_x(10:12)/h_low, h_low);
%         x_ins(10:13) = quatprod(x_ins(10:13), q_delta_omega);
% %         delta_q = [1 ; 0.5*delta_x(10:12)];
% %         x_ins(10:13) = quatprod(x_ins(10:13), delta_q);
%         x_ins(10:13) = x_ins(10:13)/norm(x_ins(10:13)); 

%         x_ins = InjectErrorState(delta_x, x_ins, g_n_hat);
        
   end
    
   
end


PlotResults;
% 
% % PLOTS
%       
% % POSITION
% figure(1)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, p_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(1,:),'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylabel('X position [m]')
% legend('True', 'Est');
% title('Position');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, p_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(2,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylabel('Y position [m]')
% legend('True', 'Est');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, p_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(3,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylabel('Z position [m]')
% legend('True', 'Est');
% grid on;
% saveas(gcf,'results/Position.jpeg')
% 
% % POSITION ERROR
% figure(2)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, p_n_nb(1,:) - ins_data(1,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-3 3]);
% ylabel('X position [m]')
% legend('Est');
% title('Position Error');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, p_n_nb(2,:) - ins_data(2,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-3 3]);
% ylabel('Y position [m]')
% legend('Error');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, p_n_nb(3,:) - ins_data(3,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-3 3]);
% ylabel('Z position [m]')
% legend('Error');
% grid on;
% saveas(gcf,'results/PositionError.jpeg')
% 
% 
% % VELOCITIES
% figure(3)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, v_n_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(4,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-20 20]);
% ylabel('X velocity [m]')
% legend('True', 'Est');
% title('Velocity');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, v_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(5,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-20 20]);
% ylabel('Y velocity [m]')
% legend('True', 'Est');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, v_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(6,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-20 20]);
% ylabel('Z velocity [m]')
% legend('True', 'Est');
% grid on;
% saveas(gcf,'results/Velocity.jpeg')
% 
% % VELOCITY ERROR
% figure(4)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, v_n_nb(1,:) - ins_data(4,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylabel('X velocity [m/s]')
% legend('Error');
% title('Velocity Error');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, v_n_nb(2,:) - ins_data(5,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylabel('Y velocity [m/s]')
% legend('Error');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, v_n_nb(3,:) - ins_data(6,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylabel('Z velocity [m/s]')
% legend('Error');
% grid on;
% saveas(gcf,'results/VelocityError.jpeg')
% 
% 
% 
% % ACCEL BIAS
% figure(5)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, bacc_b_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(7,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-20 20]);
% ylabel('X acc bias [deg]')
% legend('True', 'Est');
% title('Accelerometer bias');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, bacc_b_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(8,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-20 20]);
% ylabel('Y acc bias [deg]')
% legend('True', 'Est');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, bacc_b_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(9,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-20 20]);
% ylabel('Z acc bias [deg]')
% legend('True', 'Est');
% grid on;
% saveas(gcf,'results/AccelBiast.jpeg')
% 
% 
% % ACCEL BIAS ERROR
% figure(6)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, bacc_b_nb(1,:) - ins_data(7,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-20 20]);
% ylabel('X accel bias [m/s^2]')
% legend('Error');
% title('Accelerometer Bias Error');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, bacc_b_nb(2,:) - ins_data(8,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-20 20]);
% ylabel('Y accel bias error [m/s^2]')
% legend('Error');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, bacc_b_nb(3,:) - ins_data(9,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-20 20]);
% ylabel('Z accel bias [m/s^2]')
% legend('Error');
% grid on;
% saveas(gcf,'results/AccelBiasError.jpeg')
% 
% 
% 
% % ATTITUDE 
% figure(7)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, rad2deg*att_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
% plot(time, rad2deg*ins_data(10,:),'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-200 200]);
% ylabel('Roll angle [deg]')
% legend('True', 'Est');
% title('Attitude');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, rad2deg*att_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, rad2deg*ins_data(11,:), 'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-200 200]);
% ylabel('Pitch angle [deg]')
% legend('True', 'Est');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, rad2deg*att_n_nb(3,:),'.', 'Color', 'black', 'Linewidth', 1.5);
% plot(time, rad2deg*ins_data(12,:),'.', 'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-200 200]);
% xlabel('Time [s]');
% ylabel('Yaw angle [deg]')
% legend('True', 'Est');
% grid on;
% saveas(gcf,'results/Attitude.jpeg')
% 
% 
% % ATTITUDE ERROR
% figure(8)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, ssa(rad2deg*att_n_nb(1,:) - rad2deg*ins_data(10,:),'deg') ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-200 200]);
% ylabel('Roll angle [deg]')
% legend('Error');
% title('Attitude Error');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, ssa(rad2deg*att_n_nb(2,:) - rad2deg*ins_data(11,:), 'deg') ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-200 200]);
% ylabel('Pitch angle [deg]')
% legend('Error');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, ssa(rad2deg*att_n_nb(3,:) - rad2deg*ins_data(12,:),'deg') ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-200 200]);
% ylabel('Yaw angle [deg]')
% legend('Error');
% grid on;
% saveas(gcf,'results/AttitudeError.jpeg')
% 
% 
% 
% % GYRO BIAS
% figure(9)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, bars_b_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(13,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-1 1]);
% ylabel('Roll bias [deg/s]')
% legend('True', 'Est');
% title('Gyro bias');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, bars_b_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(14,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-1 1]);
% ylabel('Pitch bias [deg/s]')
% legend('True', 'Est');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, bars_b_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
% plot(time, ins_data(15,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
% ylim([-1 1]);
% ylabel('Yaw bias [deg/s]')
% legend('True', 'Est');
% grid on;
% saveas(gcf,'results/GyroBias.jpeg')
% 
% % GYRA BIAS ERROR
% figure(10)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, rad2deg*(bars_b_nb(1,:) - ins_data(13,:)) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-50 50]);
% ylabel('Roll bias [deg/s]')
% legend('Error');
% title('Gyro Bias Error');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, rad2deg*(bars_b_nb(2,:) - ins_data(14,:)) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-50 50]);
% ylabel('Pitch bias [deg/s]')
% legend('Error');
% grid on;
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, rad2deg*(bars_b_nb(3,:) - ins_data(15,:)) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
% ylim([-50 50]);
% ylabel('Yaw bias [deg/s]')
% legend('Error');
% grid on;
% saveas(gcf,'results/GyroBiasError.jpeg')
% 
% 
% % % Accelerometer input
% % figure(11)
% % figure(gcf);
% % subplot(3, 1, 1)
% % hold on;
% % plot(time, f_b_imu(1,:),'Color', 'black', 'Linewidth', 1.5);
% % ylabel('f_b_imu X [m/s^2]')
% % legend('Measured');
% % title('Accelerometer input');
% % grid on;
% % 
% % subplot(3, 1, 2)
% % hold on;
% % plot(time, f_b_imu(2,:),'Color', 'black', 'Linewidth', 1.5);
% % ylabel('f_b_imu y [m/s^2]')
% % legend('Measured');
% % grid on;
% % 
% % subplot(3, 1, 3)
% % hold on;
% % plot(time, f_b_imu(3,:),'Color', 'black', 'Linewidth', 1.5);
% % ylabel('f_b_imu Z [m/s^2]')
% % legend('Measured');
% % grid on;
% % % saveas(gcf,'AccelInput.jpeg')
% % 
% % 
% % 
% % % Gyrometer input
% % figure(12)
% % figure(gcf);
% % subplot(3, 1, 1)
% % hold on;
% % plot(time, rad2deg*omega_b_imu(1,:),'Color', 'black', 'Linewidth', 1.5);
% % ylabel('omega_b_imu X [m/s^2]')
% % legend('Measured');
% % title('Gyrometer input');
% % grid on;
% % 
% % subplot(3, 1, 2)
% % hold on;
% % plot(time, rad2deg*omega_b_imu(2,:),'Color', 'black', 'Linewidth', 1.5);
% % ylabel('omega_b_imu y [m/s^2]')
% % legend('Measured');
% % grid on;
% % 
% % subplot(3, 1, 3)
% % hold on;
% % plot(time, rad2deg*omega_b_imu(3,:),'Color', 'black', 'Linewidth', 1.5);
% % ylabel('omega_b_imu Z [m/s^2]')
% % legend('Measured');
% % grid on;
% % saveas(gcf,'GyroInput.jpeg')
% 
% 
% 
% % gravity error
% figure(13)
% figure(gcf);
% subplot(3, 1, 1)
% hold on;
% plot(time, g_err_data(1,:),'Color', 'black', 'Linewidth', 1.5);
% ylim([-0.00001 0.00001]);
% ylabel('Gravity X [m/s^2]')
% legend('Est');
% title('Gravity');
% grid on;
% 
% subplot(3, 1, 2)
% hold on;
% plot(time, g_err_data(2,:),'Color', 'black', 'Linewidth', 1.5);
% ylim([-0.00001 0.00001]);
% ylabel('Gravity Y [m/s^2]')
% legend('Est');
% grid on;
% 
% 
% subplot(3, 1, 3)
% hold on;
% plot(time, g_err_data(3,:),'Color', 'black', 'Linewidth', 1.5);
% ylim([9.80600 9.81000]);
% ylabel('Gravity Z [m/s^2]')
% legend('Est');
% grid on;
% 
% saveas(gcf,'results/GravityErrorState.jpeg')
% 
% % Position map
% figure(14)
% figure(gcf)
% hold on;
% plot(p_n_nb(2,:), p_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
% plot(ins_data(2,:), ins_data(1,:), 'Color', [1,165/255, 0], 'Linewidth', 1.5);
% % plot(p_n_nb(2,10*f_samp:end), p_n_nb(1,10*f_samp:end), 'Color', 'black', 'Linewidth', 1.5);
% % plot(ins_data(2,10*f_samp:end), ins_data(1,10*f_samp:end), 'Color', [1,165/255, 0], 'Linewidth', 1.5);
% xlabel('Y position [m]');
% ylabel('X position [m]');
% title('Position Plot');
% legend('True', 'Estimated');
% grid on;
% saveas(gcf,'results/PositionMap.jpeg')
% 
% 

% FUNCTIONS

%-------------------------------------------


function [x_ins_out,g_n_hat] = InjectErrorState(delta_x, x_ins, g_n_hat)
        % inject error state into nominal state
        x_ins(1:9) = x_ins(1:9) + delta_x(1:9);
        x_ins(14:16) = x_ins(14:16) + delta_x(13:15);
        g_n_hat = g_n_hat + delta_x(16:18);
%         disp(delta_x(16:18));
        h_low = 1/10;
        q_delta_omega = qbuild(delta_x(10:12)/h_low, h_low);
        x_ins(10:13) = quatprod(x_ins(10:13), q_delta_omega);
%         delta_q = [1 ; 0.5*delta_x(10:12)];
%         x_ins(10:13) = quatprod(x_ins(10:13), delta_q);
        x_ins(10:13) = x_ins(10:13)/norm(x_ins(10:13)); 
        
        x_ins_out = x_ins;
            
end

