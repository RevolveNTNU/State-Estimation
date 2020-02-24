function delta_x = ErrorStateKalman(delta_y, R_nb, Tt, init)
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