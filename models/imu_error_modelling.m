%% set up simulation parameters
clear all;
tstop = 10;
sim_dt = 1/8000; % this should be the sample time of the IMU ADC
imu_dt = 1/4000; % output data interval of the IMU
ins_dt = 1/250; % tie step that the strapdown predictor runs at
%% calculate truth trajectory
options = odeset('RelTol',1e-16);

% generate quaternions and body rates
[truth_time,truth_states]=ode113('calc_truth_deriv',[0:sim_dt:tstop],[0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0]);

% loop through the state history
body_accel = zeros(3,length(truth_time));
for index = 1:length(truth_time)
    % get the quaternion states
    quat = truth_states(index,4:7)';
    % generate a rotation matrix from earth to body
    Teb= transpose(Quat2Tbn(quat));
    % get the acceleration, subtracting gravity
    earth_accel = truth_states(index,8:10)' - [0;0;9.80665];
    % rotate into body frame
    body_accel(:,index) = Teb * earth_accel;
end

%% sample imu measurements
index_max = length(truth_time);
imu_index = 0;
delta_time_gyro = 0;
gyro_filt_state = truth_states(1,1:3)';
truth_rate_prev = truth_states(1,1:3)';
downsample_ratio = round(imu_dt/sim_dt);
imu_index_max = ceil(imu_dt/sim_dt);
imu_gyro_data = zeros(3,imu_index_max);
imu_time_data = zeros(1,imu_index_max);
imu_dt_data = zeros(1,imu_index_max);
imu_accel_data = zeros(3,imu_index_max);
% filter coefficients for gyro DLPF
fc = 250;
omegac = 2 * pi * fc * sim_dt;
K = tan(omegac / 2);
alpha = 1 + K;
a0 = 1;
a1 = -(1 - K) / alpha;
b0 = K / alpha;
b1 = K / alpha;
for index = 1:index_max;
    gyro_filt_state = (b0/a0)*truth_states(index,1:3)' + (b1/a0)*truth_rate_prev - (a1/a0)*gyro_filt_state;
    truth_rate_prev = truth_states(index,1:3)';
    % sample gyro data
    if (rem(index,downsample_ratio) == 0)
        imu_index = imu_index + 1;
        imu_gyro_data(:,imu_index) = gyro_filt_state;
        imu_time_data(imu_index) = truth_time(index);
        imu_accel_data(:,imu_index) = body_accel(:,index);
        imu_dt_data(imu_index) = imu_dt;
    end
end

%% model the IMU driver calculations
downsample_ratio = round(ins_dt/imu_dt);
method = uint8(2); % 0: addition, 1: addition + coning comp, 2: quaternion (small angle approx)
index_max = length(imu_gyro_data);
accumulated_time = 0;
ins_index = 1;
alpha = [0;0;0];
beta = [0;0;0];
last_alpha = [0;0;0];
last_delta_alpha = [0;0;0];
last_sample = [0;0;0];
gyro_rate_prev = [0;0;0];
gyro_time_prev = 0;
quat = [1;0;0;0];
delta_vel = [0;0;0];
prev_imu_accel = imu_accel_data(:,1);
for index = 1:index_max
    % integrate gyro data until the accumulated time reaches the desired
    % time step for the INS consumers
    % accumulate time
    accumulated_time = accumulated_time + imu_dt;
    % Integrate gyro data using a trapezoidal integration scheme to
    % produce the delta angle
    delta_alpha = (imu_gyro_data(:,index) + gyro_rate_prev)*imu_dt*0.5;
    gyro_rate_prev = imu_gyro_data(:,index);
    alpha = alpha + delta_alpha;
    % calculate coning corrections
    delta_beta = 0.5 * cross((last_alpha + (1.0/6.0)*last_delta_alpha),delta_alpha);
    beta = beta + delta_beta;
    last_alpha = alpha;
    % store the intermediate delta angle for use by the coning
    % correction calculation next time step
    last_delta_alpha = alpha;
    
    % try an alternative quaternion integration method
    % convert to a delta quaternion using a small angle approximation
    % and use quaternion product to accumulate rotation
    delta_quat = [1.0;0.5*delta_alpha(1);0.5*delta_alpha(2);0.5*delta_alpha(3)];
    quat = QuatMult(quat,delta_quat);
    
    % integrate accel data to get delta velocities
    delta_vel = delta_vel + (imu_accel_data(:,index) + prev_imu_accel)*imu_dt*0.5;
    prev_imu_accel = imu_accel_data(:,index);
    
    % at the prescribed downsample intervals, output the accumualted data
    % and reset the accumulator states
    if (rem(index,downsample_ratio) == 0)
        
        % output the corrected delta angle at the INS time step
        if (method == uint8(0))
            ins_delta_angle(:,ins_index) = alpha;
        elseif (method == uint8(1))
            ins_delta_angle(:,ins_index) = alpha + beta;
        elseif (method == uint8(2))
            quat=NormQuat(quat);
            ins_delta_angle(:,ins_index) = 2*[quat(2);quat(3);quat(4)];
        end
        quat = [1;0;0;0];
        
        % output the delta velocity at the INS time step
        ins_delta_vel(:,ins_index) = delta_vel;
        
        % output the time stamp of the INS delta angles - use the time from
        % the most recent gyro measurement
        ins_time(ins_index) = imu_time_data(index);
        ins_time_delta(ins_index) = accumulated_time;
        
        % reset the accumulated values
        alpha = [0;0;0];
        beta = [0;0;0];
        delta_vel = [0;0;0];
        accumulated_time = 0;
        
        % increment the storage index
        ins_index = ins_index + 1;
    end
end

%% model the strapdown calculations
index_max = length(ins_delta_angle);
quat = [1;0;0;0];
vel = [0;0;0];
vel_prev = vel;
pos = [0;0;0];
Tbn_prev = Quat2Tbn(quat);
for index = 1:index_max
    % Convert the rotation vector to its equivalent quaternion
    deltaQuat = RotToQuat(ins_delta_angle(:,index));
    
    % Update the quaternions by rotating from the previous attitude through
    % the delta angle rotation quaternion
    quat = QuatMult(quat,deltaQuat);
    quat = NormQuat(quat);
    
    % calculate the rotation matrix from Body to earth frame
    Tbn = Quat2Tbn(quat);
    
    % rotate the delta velocity data into earth frame and correct for
    % gravity
    %del_vel_earth = 0.5*(Tbn * ins_delta_vel(:,index) + Tbn_prev * ins_delta_vel(:,index)) + [0;0;9.80665]*ins_dt;
    del_vel_earth = Tbn * ins_delta_vel(:,index) + [0;0;9.80665]*ins_time_delta(index);
    Tbn_prev = Tbn;
    
    % sum to give velocity estimate
    vel = vel + del_vel_earth;
    
    % integrate to give position
    pos = pos + (vel + vel_prev)*0.5*ins_time_delta(index);
    vel_prev = vel;
    
    % convert to Euler angles for visualisation
    ins_euler(:,index) = QuatToEul(quat);
    ins_vel(:,index) = vel;
    ins_pos(:,index) = pos;
    
end

%% get truth data at INS calculation times
for index = 1:index_max
    ref_index = find(truth_time == ins_time(index));
    quat = truth_states(ref_index,4:7)';
    ref_euler(:,index) = QuatToEul(quat);
    ref_time(index) = ins_time(index);
    ref_vel(:,index) = truth_states(ref_index,11:13);
    ref_pos(:,index) = truth_states(ref_index,14:16);
end

%% plot the data
figure;
subplot(3,1,1);
plot(ins_time,180/pi*ins_euler(1,:));
title('Euler angle comparison');
xlabel('time(sec)');ylabel('roll (deg)');grid on;
hold on;
plot(ref_time,180/pi*ref_euler(1,:));
hold off;
subplot(3,1,2);
plot(ins_time,180/pi*ins_euler(2,:));
xlabel('time(sec)');ylabel('pitch (deg)');grid on;
hold on;
plot(ref_time,180/pi*ref_euler(2,:));
hold off;
subplot(3,1,3);
plot(ins_time,180/pi*ins_euler(3,:));
xlabel('time(sec)');ylabel('yaw (deg)');grid on;
hold on;
plot(ref_time,180/pi*ref_euler(3,:));
hold off;

figure;
subplot(3,1,1);
plot(ins_time,ins_vel(1,:));
title('velocity comparison');
xlabel('time(sec)');ylabel('vx (m/s)');grid on;
hold on;
plot(ref_time,ref_vel(1,:));
hold off;
subplot(3,1,2);
plot(ins_time,ins_vel(2,:));
xlabel('time(sec)');ylabel('vy (m/s)');grid on;
hold on;
plot(ref_time,ref_vel(2,:));
hold off;
subplot(3,1,3);
plot(ins_time,ins_vel(3,:));
xlabel('time(sec)');ylabel('vz (m/s)');grid on;
hold on;
plot(ref_time,ref_vel(3,:));
hold off;