%% Derive an EKF to estimate the NED position of the origin of a range beacon system using range to anchor and laser range to deck measurements
clear all;
close all;

syms bcn_pos_n bcn_pos_e bcn_pos_d real; % NED position of beacons relative to landing deck origin in earth frame 
syms veh_pos_n veh_pos_e veh_pos_d real; % NED position of vehicle in earth frame
syms deck_pos_n deck_pos_e deck_pos_d real; % NED position of deck origin in earth frame
syms deck_vel_n deck_vel_e deck_vel_d real; % NED velocity of deck origin in earth frame
syms dt real;
syms R_RNG real;

stateVector =    [deck_vel_n                   ; deck_vel_e                   ; deck_vel_d ;...
                  deck_pos_n                   ; deck_pos_e                   ; deck_pos_d];
newStateVector = [deck_vel_n                   ; deck_vel_e                   ; deck_vel_d ;...
                  deck_pos_n + deck_vel_n * dt ; deck_pos_e + deck_vel_e * dt ; veh_pos_d + deck_vel_d * dt];

% derive the state transition matrix
F = jacobian(newStateVector, stateVector);

% define a symbolic covariance matrix using strings to represent
% '_l_' to represent '( '
% '_c_' to represent ,
% '_r_' to represent ')'
% these can be substituted later to create executable code
for rowIndex = 1:6
    for colIndex = 1:6
        eval(['syms OP_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['P(',num2str(rowIndex),',',num2str(colIndex), ') = OP_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

% Derive the predicted covariance matrix using the standard equation
PP = F*P*transpose(F);

%% derive equations for fusion of range to anchor measurements
rng_pred = sqrt((veh_pos_n - deck_pos_n - bcn_pos_n)^2 + (veh_pos_e - deck_pos_e - bcn_pos_e)^2 + (veh_pos_d - deck_pos_d - bcn_pos_d)^2);% derive equations for a simple 3 state estimator

H_RNG = jacobian(rng_pred,stateVector); % measurement Jacobian
[H_RNG,SH_RNG]=OptimiseAlgebra(H_RNG,'SH_RNG');

% calculate the Kalman gain matrix and optimise algebra
K_RNG = (P*transpose(H_RNG))/(H_RNG*P*transpose(H_RNG) + R_RNG);
[K_RNG,SK_RNG]=OptimiseAlgebra(K_RNG,'SK_RNG');

%% derive fusion of Z axis laser range finder
syms qd0 qd1 qd2 qd3 real % quaternions defining attitude of platform XYZ axes relative to local NED
syms qv0 qv1 qv2 qv3 real % quaternions defining attitude of vehicle XYZ axes relative to local NED

% derive the deck to nav direction cosine matrix
Tdn = Quat2Tbn([qd0,qd1,qd2,qd3]);

% derive the vehicle to nav direction cosine matrix
Tvn = Quat2Tbn([qv0,qv1,qv2,qv3]);

% define vehicle position in nav reference frame
pos_veh_nav = [veh_pos_n ; veh_pos_e ; veh_pos_d];

% rotate into deck frame
pos_veh_deck = transpose(Tdn) * pos_veh_nav;

% define deck surface in nav frame using position of deck origin and deck normal (Z) unit vector
deck_normal_vec = Tdn * [0;0;1];
deck_origin = [deck_pos_n;deck_pos_e;deck_pos_d];

% define range finder unit vector
veh_rngfnd_vec = Tvn * [0;0;1];

% define position along range finder vector as a function of range
% measurement
syms p_n p_e p_d real % NED position of intersection of laser and deck
syms rng_laser real % path length along range finder vector
p_n = veh_pos_n + veh_rngfnd_vec(1) * rng_laser;
p_e = veh_pos_e + veh_rngfnd_vec(2) * rng_laser;
p_d = veh_pos_d + veh_rngfnd_vec(3) * rng_laser;

% find the range where the laser beam and deck plane intersect
rng_laser = solve(deck_normal_vec(1)*(p_n - deck_origin(1)) + ...
      deck_normal_vec(2)*(p_e - deck_origin(2)) + ...
      deck_normal_vec(3)*(p_d - deck_origin(3)),'rng_laser');

H_RNG2 = jacobian(rng_laser,stateVector); % measurement Jacobian
[H_RNG2,SH_RNG2]=OptimiseAlgebra(H_RNG2,'SH_RNG2');

% calculate the Kalman gain matrix and optimise algebra
syms R_LASER real % laser range finder observation variance
K_RNG2 = (P*transpose(H_RNG2))/(H_RNG2*P*transpose(H_RNG2) + R_LASER);
[K_RNG2,SK_RNG2]=OptimiseAlgebra(K_RNG2,'SK_RNG2');

%% Save output and convert to m and c code fragments
nStates = length(PP);
fileName = strcat('SymbolicOutput',int2str(nStates),'.mat');
save(fileName);
SaveScriptCode(nStates);
ConvertToM(nStates);
ConvertToC(nStates);
