function [state_deriv] = calc_truth_deriv(time,state)

% state vector = [P;Q;R;q0;q1;q2;q3;q4]

% model a linearly increasing harmonic angular acceleration about X and Y
% where X and Y are 90 degrees out of phase
freq = 50*2*pi;
amplitude = 10*time;
ang_accel = [amplitude*sin(freq*time);amplitude*cos(freq*time);0];

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

% assemble the state derivative vector 
state_deriv = [ang_accel;quat_deriv];

end

