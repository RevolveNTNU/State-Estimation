function data = LoadVCUData(folder_name,start_index, end_index)

data = ImportAnalyzeData(folder_name);
data.vcu_INS_ax = data.vcu_INS_ax(start_index:end_index);
data.vcu_INS_ay = data.vcu_INS_ay(start_index:end_index);
data.vcu_INS_az = data.vcu_INS_az(start_index:end_index);

coeff = ones(1,20)/20;
data.vcu_INS_ax = filter(coeff, 1, data.vcu_INS_ax);
data.vcu_INS_ay = filter(coeff, 1, data.vcu_INS_ay);
data.vcu_INS_az = filter(coeff, 1, data.vcu_INS_az);

data.vcu_INS_vx = data.vcu_INS_vx(start_index:end_index);
data.vcu_INS_vy = data.vcu_INS_vy(start_index:end_index);
data.vcu_INS_vz = data.vcu_INS_vz(start_index:end_index);

data.vcu_INS_yaw = data.vcu_INS_yaw(start_index:end_index);
data.vcu_INS_yaw_rate = data.vcu_INS_yaw_rate(start_index:end_index);
data.vcu_INS_roll = data.vcu_INS_roll(start_index:end_index);
data.vcu_INS_roll_rate = data.vcu_INS_roll_rate(start_index:end_index);
data.vcu_INS_pitch = data.vcu_INS_pitch(start_index:end_index);
data.vcu_INS_pitch_rate = data.vcu_INS_pitch_rate(start_index:end_index);

end

