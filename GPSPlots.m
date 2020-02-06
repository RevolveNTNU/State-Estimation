

l0 = SKID.vcu_GNSS_longitude(10);
mu0 = SKID.vcu_GNSS_latitude(10);
h_ref = SKID.vcu_GNSS_altitude(10);

gps_f = 5; %hz
gps_h = 1/gps_f;
gps_simtime = 30;
gps_N = gps_simtime/gps_h;

gps_time_data = zeros(1,gps_N);
gps_pos_data = zeros(3,gps_N);

count = 20;
index = 1;
for i = 2:22800 
   count = count + 1;
   if (count >= 20)
       count = 0;
       gps_t = (i-1) * gps_h;
       gps_time_data(i) = gps_t;

       l = SKID.vcu_GNSS_longitude(i);
       mu = SKID.vcu_GNSS_latitude(i);
       gps_h = SKID.vcu_GNSS_altitude(i);

       [x,y,z] = llh2flat(l,mu,gps_h,l0,mu0,h_ref);
       gps_pos_data(1,index) = 0.2 * -x;
       gps_pos_data(2,index) = 0.2 * -y;
       gps_pos_data(3,index) = 0.2 * -z;
       
       index = index+1;
   end
end

%PLOTS

gps_t = gps_time_data;
gps_xpos = gps_pos_data(1,:);
gps_ypos = gps_pos_data(2,:);
gps_zpos = gps_pos_data(3,:);

figure(10)
figure(gcf);
subplot(1,1,1)
hold on;
plot(gps_xpos, gps_ypos, 'Color', 'black', 'Linewidth', 1.5);
xlabel('X position [m]')
ylabel('Y position [m]')
legend('GPS position');


figure(11)
figure(gcf);
subplot(1,1,1)
hold on;
plot(SKID.vcu_GNSS_longitude, SKID.vcu_GNSS_latitude, 'Color', 'black', 'Linewidth', 1.5);
xlabel('X position [m]')
ylabel('Y position [m]')
legend('GPS position');
% 
% subplot(3,1,2)
% hold on;
% plot(gps_t, gps_ypos, 'Color', 'black', 'Linewidth', 1.5);
% ylabel('Y position [m]')
% legend('GPS position y');
% 
% subplot(3,1,3)
% hold on;
% plot(gps_t, gps_zpos, 'Color', 'black', 'Linewidth', 1.5);
% ylabel('Z position [m]')
% legend('GPS position z');
% 
