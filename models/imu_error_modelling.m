%% set up simulation parameters
clear all;
tstop = 10;
sim_dt = 1e-4;
imu_dt = 1e-3;
ins_dt = 4e-3;
%% calculate truth trajectory
options = odeset('RelTol',1e-12);
[truth_time,truth_states]=ode113('calc_truth_deriv',[0:sim_dt:tstop],[0;0;0;1;0;0;0]);

%% sample gyro rate measurements
index_max = length(truth_time);
imu_index = 0;
delta_time_gyro = 0;
for index = 1:index_max;
    % sample gyro data
    delta_time_gyro = delta_time_gyro + sim_dt;
    if (delta_time_gyro >= imu_dt)
        imu_index = imu_index + 1;
        gyro_rate_data(:,imu_index) = truth_states(index,1:3)';
        gyro_time_data(imu_index) = truth_time(index);
        delta_time_gyro = 0;
    end
end

%% model the IMU driver calculations
index_max = length(gyro_rate_data);
accumulated_time = 0;
ins_index = 1;
alpha = [0;0;0];
beta = [0;0;0];
last_alpha = [0;0;0];
last_delta_alpha = [0;0;0];
last_sample = [0;0;0];
gyro_rate_prev = [0;0;0];
gyro_time_prev = 0;
for index = 1:index_max
    % integrate gyro data until the accumulated time reaches the desired
    % time step for the INS consumers
    if (accumulated_time < ins_dt)
        % accumulate time
        delta_time = gyro_time_data(index) - gyro_time_prev;
        gyro_time_prev = gyro_time_data(index);
        accumulated_time = accumulated_time + delta_time;
        % Integrate gyro data using a trapezoidal integration scheme to
        % produce the delta angle
        delta_alpha = (gyro_rate_data(:,index) + gyro_rate_prev)*delta_time*0.5;
        gyro_rate_prev = gyro_rate_data(:,index);
        alpha = alpha + delta_alpha;
        % calculate coning corrections
        delta_beta = 0.5 * cross((last_alpha + (1.0/6.0)*last_delta_alpha),delta_alpha);
        beta = beta + delta_beta;
        last_alpha = alpha;
        % store the intermediate delta angle for use by the coning
        % correction calculation next time step
        last_delta_alpha = alpha;
    else
        % output the corrected delta angle at the INS time step
        ins_delta_angle(:,ins_index) = alpha + beta;
        % output the time stamp of the INS delta angles - use the time from
        % the most recent gyro measurement
        ins_time(ins_index) = gyro_time_data(index);
        % reset the accumulated values
        alpha = [0;0;0];
        beta = [0;0;0];
        accumulated_time = 0;
        % increment the storage index
        ins_index = ins_index + 1;
    end
end

%% model the strapdown calculations
index_max = length(ins_delta_angle);
quat = [1;0;0;0];
for index = 1:index_max
    % Convert the rotation vector to its equivalent quaternion
    deltaQuat = RotToQuat(ins_delta_angle(:,index));
    
    % Update the quaternions by rotating from the previous attitude through
    % the delta angle rotation quaternion
    quat = QuatMult(quat,deltaQuat);
    quat = NormQuat(quat);
    
    % convert to Euler angles for visualisation
    ins_euler(:,index) = QuatToEul(quat);
end

%% get truth data at INS calculation times
for index = 1:index_max
    ref_index = find(truth_time == ins_time(index));
    quat = truth_states(ref_index,4:7)';
    ref_euler(:,index) = QuatToEul(quat);
    ref_time(index) = ins_time(index);
end

%% plot the data
figure;
subplot(3,1,1);
plot(ins_time,180/pi*ins_euler(1,:));
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