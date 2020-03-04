function delta_q = qbuild(omega,h)
%QBUILD Summary of this function goes here
%   Detailed explanation goes here

alpha_x = omega(1);
alpha_y = omega(2);
alpha_z = omega(3);

unit_vec_x = [1 0 0]';
unit_vec_y = [0 1 0]';
unit_vec_z = [0 0 1]';


qx = [cos(0.5 * h * alpha_x) ; sin(0.5 * h * alpha_x) * unit_vec_x];
qy = [cos(0.5 * h * alpha_y) ; sin(0.5 * h * alpha_y) * unit_vec_y];
qz = [cos(0.5 * h * alpha_z) ; sin(0.5 * h * alpha_z) * unit_vec_z];

q = quatprod(qx,qy);
q = quatprod(q,qz);

delta_theta = omega*h;
delta_theta_norm = norm(delta_theta);

if(delta_theta_norm > 1e-10) 
    delta_q = [ cos(delta_theta_norm * 0.5) ; sin(delta_theta_norm * 0.5) * (delta_theta/delta_theta_norm)];
else
    delta_q = [1 0 0 0]';
end

end

