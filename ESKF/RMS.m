function [rms_tot, rms_lane, rms_skid] = RMS(true, est)

sum_tot = 0;
sum_lane = 0;
sum_skid = 0;

for i = 2000:length(true)
    sum_tot = sum_tot + (true(1:3,i) - est(1:3,i)).^2;
end

for i = 7701:8250
    sum_lane = sum_lane + (true(1:3,i) - est(1:3,i)).^2;
end

for i = 8251:length(true)
    sum_skid = sum_skid + (true(1:3,i) - est(1:3,i)).^2;
end

rms_tot = sum_tot / length(true);
rms_lane = sum_lane / (8250 - 7701);
rms_skid = sum_skid / (length(true) - 8251);


end

