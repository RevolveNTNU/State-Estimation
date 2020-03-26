clc
clear all

%%

z31 = zeros(3,1);
z_151 = zeros(15,1);

%% pos
syms px_n py_n pz_n real
p_hat = [px_n; py_n; pz_n];
syms d_px_n d_py_n d_pz_n real
d_p = [d_px_n; d_py_n; d_pz_n];
p = p_hat + d_p

%% vel 

syms vx_n vy_n vz_n real
v_hat = [vx_n; vy_n; vz_n];
syms d_vx_n d_vy_n d_vz_n real
d_v = [d_vx_n; d_vy_n; d_vz_n];

v = v_hat + d_v

%% acc
syms b_accx_b b_accy_b b_accz_b real
b_acc_b_hat = [b_accx_b; b_accy_b; b_accz_b];
syms d_b_accx_b d_b_accy_b d_b_accz_b real
d_b_acc = [d_b_accx_b; d_b_accy_b; d_b_accz_b];

b_acc_b = b_acc_b_hat + d_b_acc


%% att
syms d_theta_x d_theta_y d_theta_z real
d_theta = [d_theta_x; d_theta_y; d_theta_z];
d_R = eye(3) + Smtrx(d_theta)
%syms q_w q_x q_y q_z real
%q = [q_w; q_x; q_y; q_z];
%R_hat = quat2rotMat_fast(q )
syms R11 R12 R13 R21 R22 R23 R31 R32 R33 real
R_hat = [...
    R11 R12 R13; 
    R21 R22 R23;
    R31 R32 R33 ]
R = R_hat*d_R;

%% ars
syms b_arsx_b b_arsy_b b_arsz_b real
b_ars_b_hat = [b_arsx_b; b_arsy_b; b_arsz_b];
syms d_b_arsx_b d_b_arsy_b d_b_arsz_b real
d_b_ars = [d_b_arsx_b; d_b_arsy_b; d_b_arsz_b];

b_ars_b = b_ars_b_hat + d_b_ars

%% error state def
d_x = [d_p; d_v; d_b_acc; d_theta; d_b_ars];

%% imu
% syms ox_b oy_b oz_b real
% omega_nb = [ox_b oy_b oz_b];
% syms wx_ars_b wy_ars_b wz_ars_b real
% w_ars_b = [wx_ars_b; wy_ars_b; wz_ars_b]
syms ox_imu_b oy_imu_b oz_imu_b real
omega_imu_b = [ox_imu_b; oy_imu_b; oz_imu_b]

%% lever arm
syms rx_b ry_b rz_b real
r_b = [rx_b; ry_b; rz_b ]

%% simple measurement
H_alloc = [1 0 0];
y_gss = H_alloc*( R'*v + Smtrx( omega_imu_b - b_ars_b)*r_b );
y_hat = subs(y_gss, d_x, z_151 )
H = jacobian( y_gss, d_x );
H = simplify( subs( H, d_x, zeros(15,1) ));
H_pos = H(:,1:3)
H_vel = H(:,4:6)
H_acc = H(:,7:9)
H_att = H(:,10:12)
H_ars = H(:,13:15)

% matriseregne sjekk, svar skal bli null
delta_H_pos = H(:,1:3) - H_alloc*z31
delta_H_vel = H(:,4:6) - H_alloc*R_hat'
delta_H_acc = H(:,7:9) - H_alloc*z31 
delta_H_att = H(:,10:12) - H_alloc*Smtrx( R_hat'*v_hat )
delta_H_ars = H(:,13:15) - [1 0 0]*Smtrx(r_b)

%% measurement
H_alloc = [1 0 0; 0 1 0; 0 0 0];
velocity_gss = H_alloc*( R'*v + Smtrx( omega_imu_b - b_ars_b)*r_b );
y_gss = norm( velocity_gss );
y_hat = simplify( subs(y_gss, d_x, z_151 ) );
H = jacobian( y_gss, d_x );
H = simplify( subs( H, d_x, zeros(15,1) ) )

%
H_pos = H(:,1:3)
H_vel = H(:,4:6)
H_acc = H(:,7:9)
H_att = H(:,10:12)
H_ars = H(:,13:15)

% se på forskjellen til 
H_alloc = [1 0 0];
y_gss = H_alloc*( R'*v + Smtrx( omega_imu_b - b_ars_b)*r_b ) 
% og 
H_alloc = [1 0 0; 0 1 0; 0 0 0];
velocity_gss = H_alloc*( R'*v + Smtrx( omega_imu_b - b_ars_b)*r_b );
y_gss = norm( velocity_gss );
% i simulatoren