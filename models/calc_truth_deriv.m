function [state_deriv] = calc_truth_deriv(time,state)

% state vector = [P;Q;R;q0;q1;q2;q3;ax;ay;az;vx;vy;vz;px;py;pz]

% model a linearly increasing harmonic angular acceleration about X and Y
% where X and Y are 90 degrees out of phase
freq = 250*2*pi;
amplitude = 10*time;
% input a yaw acceleration doublet
if (time <= 1.0) 
    ang_acc_z = 10.0*time;
elseif (time <= 3.0)
    ang_acc_z = 20.0 - 10.0*time;
elseif (time <= 4.0)
    ang_acc_z = - 40.0 + 10.0*time;
else
    ang_acc_z = 0.0;
end
ang_accel = [amplitude*sin(freq*time);amplitude*cos(freq*time);ang_acc_z];

% calculate the transfer matrix from quaternions to quaternion rates
omega_mat = [ ...
    0,state(1),state(2),state(3); ...
    -state(1),0,-state(3),state(2); ...
    -state(2),state(3),0,-state(1); ...
    -state(3),-state(2),state(1),0  ...
    ];

% quaternion vector
quat = [state(4);state(5);state(6);state(7)];

% quaternion normalisation error
epsilon = 1 - sum(quat.*quat);

% calculate the quaternion derivative with a correction for normalisation
% errors
quat_deriv = - 0.5 * omega_mat * quat + quat * epsilon * 1e-6;

% calculate the acceleration derivative harmonic vibration
accel_freq_x = 135*2*pi;
accel_freq_y = 115*2*pi;
accel_freq_z = 126*2*pi;
acc_deriv = 10*[accel_freq_x*sin(accel_freq_x*time);accel_freq_y*sin(accel_freq_y*time);accel_freq_z*sin(accel_freq_z*time)] - 10.0*[state(11);state(12);state(13)] - 10.0*[state(8);state(9);state(10)];

% calculate the position derivative
vel_deriv = [state(8);state(9);state(10)];

% calculate the position derivative
pos_deriv = [state(11);state(12);state(13)];

% assemble the state derivative vector 
state_deriv = [ang_accel;quat_deriv;acc_deriv;vel_deriv;pos_deriv];

end

