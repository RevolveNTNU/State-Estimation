% PLOTS

zero_ref = zeros(1,length(p_n_nb));

% POSITION
figure(1)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, p_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(1,:),'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylabel('X position [m]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
title('Position');
grid minor;

subplot(3, 1, 2)
hold on;
plot(time, p_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(2,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylabel('Y position [m]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
grid minor;

subplot(3, 1, 3)
hold on;
plot(time, p_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(3,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylabel('Z position [m]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
grid minor;
saveas(gcf,'results_no-vec-no-gnss1/Position.jpeg')

% POSITION ERROR
figure(2)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, p_n_nb(1,:) - ins_data(1,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-1 1]);
ylabel('X position [m]')
xlabel('Time [s]')
legend('Zero ref.', 'Error','Orientation','horizontal');
title('Position Error');
grid minor;

subplot(3, 1, 2)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, p_n_nb(2,:) - ins_data(2,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);ylim([-1 1]);
ylabel('Y position [m]')
xlabel('Time [s]')
legend('Zero ref.', 'Error','Orientation','horizontal');
grid minor;

subplot(3, 1, 3)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, p_n_nb(3,:) - ins_data(3,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-1 1]);
ylabel('Z position [m]')
xlabel('Time [s]')
legend('Zero ref.', 'Error','Orientation','horizontal');
grid minor;
saveas(gcf,'results_no-vec-no-gnss1/PositionError.jpeg')


% VELOCITIES
figure(3)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, v_n_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(4,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-15 15]);
ylabel('X velocity [m]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
title('Velocity');
grid minor;

subplot(3, 1, 2)
hold on;
plot(time, v_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(5,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-15 15]);
ylabel('Y velocity [m]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
grid minor;

subplot(3, 1, 3)
hold on;
plot(time, v_n_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(6,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-15 15]);
ylabel('Z velocity [m]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
grid minor;
saveas(gcf,'results_no-vec-no-gnss1/Velocity.jpeg')

% VELOCITY ERROR
figure(4)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, v_n_nb(1,:) - ins_data(4,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-1 1]);
ylabel('X velocity [m/s]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
title('Velocity Error');
grid minor;

subplot(3, 1, 2)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, v_n_nb(2,:) - ins_data(5,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-1 1]);
ylabel('Y velocity [m/s]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
grid minor;

subplot(3, 1, 3)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, v_n_nb(3,:) - ins_data(6,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-1 1]);
ylabel('Z velocity [m/s]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
grid minor;
saveas(gcf,'results_no-vec-no-gnss1/VelocityError.jpeg')



% ACCEL BIAS
figure(5)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, bacc_b_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(7,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-1 1]);
ylabel('X acc bias [deg]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
title('Accelerometer bias');
grid minor;

subplot(3, 1, 2)
hold on;
plot(time, bacc_b_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(8,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-1 1]);
ylabel('Y acc bias [deg]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
grid minor;

subplot(3, 1, 3)
hold on;
plot(time, bacc_b_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(9,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-1 1]);
ylabel('Z acc bias [deg]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
grid minor;
saveas(gcf,'results_no-vec-no-gnss1/AccelBiast.jpeg')


% ACCEL BIAS ERROR
figure(6)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, bacc_b_nb(1,:) - ins_data(7,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-1 1]);
ylabel('X accel bias [m/s^2]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
title('Accelerometer Bias Error');
grid minor;

subplot(3, 1, 2)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, bacc_b_nb(2,:) - ins_data(8,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-1 1]);
ylabel('Y accel bias error [m/s^2]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
grid minor;

subplot(3, 1, 3)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, bacc_b_nb(3,:) - ins_data(9,:) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-1 1]);
ylabel('Z accel bias [m/s^2]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
grid minor;
saveas(gcf,'results_no-vec-no-gnss1/AccelBiasError.jpeg')



% ATTITUDE 
figure(7)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, rad2deg*att_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
plot(time, rad2deg*ins_data(10,:),'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-20 20]);
ylabel('Roll angle [deg]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
title('Attitude');
grid minor;

subplot(3, 1, 2)
hold on;
plot(time, rad2deg*att_n_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, rad2deg*ins_data(11,:), 'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-20 20]);
ylabel('Pitch angle [deg]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
grid minor;

subplot(3, 1, 3)
hold on;
plot(time, rad2deg*att_n_nb(3,:),'.', 'Color', 'black', 'Linewidth', 1.5);
plot(time, rad2deg*ins_data(12,:),'.', 'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-200 200]);
xlabel('Time [s]');
ylabel('Yaw angle [deg]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
grid minor;
saveas(gcf,'results_no-vec-no-gnss1/Attitude.jpeg')


% ATTITUDE ERROR
figure(8)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, ssa(rad2deg*att_n_nb(1,:) - rad2deg*ins_data(10,:),'deg') ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-5 5]);
ylabel('Roll angle [deg]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
title('Attitude Error');
grid minor;

subplot(3, 1, 2)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, ssa(rad2deg*att_n_nb(2,:) - rad2deg*ins_data(11,:), 'deg') ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-5 5]);
ylabel('Pitch angle [deg]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
grid minor;

subplot(3, 1, 3)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, ssa(rad2deg*att_n_nb(3,:) - rad2deg*ins_data(12,:),'deg') ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-5 5]);
ylabel('Yaw angle [deg]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
grid minor;
saveas(gcf,'results_no-vec-no-gnss1/AttitudeError.jpeg')



% GYRO BIAS
figure(9)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, bars_b_nb(1,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(13,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-0.1 0.1]);
ylabel('Roll bias [deg/s]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
title('Gyro bias');
grid minor;

subplot(3, 1, 2)
hold on;
plot(time, bars_b_nb(2,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(14,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-0.1 0.1]);
ylabel('Pitch bias [deg/s]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
grid minor;

subplot(3, 1, 3)
hold on;
plot(time, bars_b_nb(3,:),'Color', 'black', 'Linewidth', 1.5);
plot(time, ins_data(15,:) ,'Color', [1,165/255, 0], 'Linewidth', 1.5);
ylim([-0.1 0.1]);
ylabel('Yaw bias [deg/s]')
xlabel('Time [s]')
legend('True', 'Est', 'Orientation','horizontal');
grid minor;
saveas(gcf,'results_no-vec-no-gnss1/GyroBias.jpeg')

% GYRA BIAS ERROR
figure(10)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, rad2deg*(bars_b_nb(1,:) - ins_data(13,:)) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-1 1]);
ylabel('Roll bias [deg/s]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
title('Gyro Bias Error');
grid minor;

subplot(3, 1, 2)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, rad2deg*(bars_b_nb(2,:) - ins_data(14,:)) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-1 1]);
ylabel('Pitch bias [deg/s]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
grid minor;

subplot(3, 1, 3)
hold on;
plot(time, zero_ref, '--', 'Color', 'black', 'Linewidth', 0.5);
plot(time, rad2deg*(bars_b_nb(3,:) - ins_data(15,:)) ,'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 0.5);
ylim([-1 1]);
ylabel('Yaw bias [deg/s]')
xlabel('Time [s]')
legend('Zero ref.', 'Error', 'Orientation','horizontal');
grid minor;
saveas(gcf,'results_no-vec-no-gnss1/GyroBiasError.jpeg')


% Gravity (estimate)
figure(13)
figure(gcf);
subplot(3, 1, 1)
hold on;
plot(time, g_err_data(1,:),'Color', 'black', 'Linewidth', 1.5);
% ylim([-0.00001 0.00001]);
ylabel('Gravity X [m/s^2]')
xlabel('Time [s]')
legend('Est');
title('Gravity');
grid minor;

subplot(3, 1, 2)
hold on;
plot(time, g_err_data(2,:),'Color', 'black', 'Linewidth', 1.5);
% ylim([-0.00001 0.00001]);
ylabel('Gravity Y [m/s^2]')
xlabel('Time [s]')
legend('Est');
grid minor;


subplot(3, 1, 3)
hold on;
plot(time, g_err_data(3,:),'Color', 'black', 'Linewidth', 1.5);
% ylim([9.80600 9.81000]);
ylabel('Gravity Z [m/s^2]')
xlabel('Time [s]')
legend('Est');
grid minor;

saveas(gcf,'results_no-vec-no-gnss1/GravityErrorState.jpeg')

% Position map
figure(14)
figure(gcf)
hold on;
plot(p_n_nb(2,:), p_n_nb(1,:), 'Color', 'black', 'Linewidth', 1.5);
plot(ins_data(2,:), ins_data(1,:), 'Color', [1,165/255, 0], 'Linewidth', 1.5);
% plot(p_n_nb(2,10*f_samp:end), p_n_nb(1,10*f_samp:end), 'Color', 'black', 'Linewidth', 1.5);
% plot(ins_data(2,10*f_samp:end), ins_data(1,10*f_samp:end), 'Color', [1,165/255, 0], 'Linewidth', 1.5);
xlabel('Y position [m]');
ylabel('X position [m]');
title('Position Plot');
legend('True', 'Est', 'Orientation','horizontal');
grid minor;
saveas(gcf,'results_no-vec-no-gnss1/PositionMap.jpeg')
