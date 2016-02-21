% IMPORTANT - This script requires the Matlab symbolic toolbox
% Derivation of EKF equations for estimation of terrain height offset
% Author:  Paul Riseborough

% State vector:

% vehicle vertical velocity in local NED frame
% vehicle vertical position in loal NED frame
% terrain vertical position in local NED frame


% Observations:

% vehicle vertical position estimate from the nav system relative to local NED frame
% Range measurement aligned with the Z body axis assuming a flat earth model

% Control inputs parameters:

% vertical rate of change of velocity in local NED frame

clear all;

%% define symbolic variables and constants
% states
syms VD real % down velocity : m/sec
syms PD real % down position : m
syms TPD real % down position of terrain : m
% state noise
syms R_TPD real % process noise variance of terrain position ; m^2
% observations
syms meaHAGL real; % height above ground measurement : m
syms meaPD real; % vertical position measurement : m
syms R_HAGL real % variance of range finder measurement : m^2
syms R_PD real; % variance of vertical position measurement : m^2
% control inputs
syms AD real % downwards rate of change of velocity : m/sec^2
syms R_AD real % noise variance of AD : (m/sec^2)^2
% parameters
syms dt real % time step : sec
nStates = 3;

%% define the process equations

% define the state vector
stateVector = [VD;PD;TPD];

% define the process equations
newVD = VD + dt*AD;
newPD = PD + dt*VD;
newTPD = TPD;
processEqns = [newVD;newPD;newTPD];

%% derive the state transition matrix

% derive the state transition matrix
F = jacobian(processEqns, stateVector);

%% derive the covariance prediction equation
% This reduces the number of floating point operations by a factor of 6 or
% more compared to using the standard matrix operations in code

% Define the control (disturbance) vector. Error growth in the vehicle
% vertical velocity and position estimates are assumed to be driven by 
% 'noise' in the vertical acceleration.
distVector = [AD];

% derive the control(disturbance) influence matrix
G = jacobian(processEqns, distVector);

% derive the state error matrix for control inputs
imuNoise = diag([R_AD]);
Q = G*imuNoise*transpose(G);

% define a symbolic covariance matrix using strings to represent 
% '_lp_' to represent '( '
% '_c_' to represent ,
% '_rp_' to represent ')' 
% these can be substituted later to create executable code
for rowIndex = 1:nStates
    for colIndex = 1:nStates
        eval(['syms OP_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['P(',num2str(rowIndex),',',num2str(colIndex), ') = OP_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

% Derive the predicted covariance matrix using the standard equation
PP = F*P*transpose(F) + Q + diag([0;0;R_TPD]);

% Collect common expressions to optimise processing
[PP,SPP]=OptimiseAlgebra(PP,'SPP');

%% derive equations for fusion of height above ground
meaHAGL = TPD - PD;
H_HAGL = jacobian(meaHAGL,stateVector); % measurement Jacobian

% calculate the Kalman gain matrix and optimise algebra
K_HAGL = (P*transpose(H_HAGL))/(H_HAGL*P*transpose(H_HAGL) + R_HAGL);
[K_HAGL,SK_HAGL]=OptimiseAlgebra(K_HAGL,'SK_HAGL');

%% derive equations for fusion of vehice position measurements
meaPD = PD;
H_PD = jacobian(meaPD,stateVector); % measurement Jacobian

% calculate the Kalman gain matrix and optimise algebra
K_PD = (P*transpose(H_PD))/(H_PD*P*transpose(H_PD) + R_PD);
[K_PD,SK_PD]=OptimiseAlgebra(K_PD,'SK_PD');

%% Save output
fileName = strcat('SymbolicOutput3.mat');
save(fileName);

%% Write equations to text and convert to m and c code fragments
SaveScriptCodeTerrainEstimator3;
ConvertToM(3);
ConvertToC(3);
