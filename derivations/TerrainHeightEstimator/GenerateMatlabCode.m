% IMPORTANT - This script requires the Matlab symbolic toolbox
% Derivation of EKF for estimation of terrain height offset
% Author:  Paul Riseborough
% Last Modified: 4 Jan 2015

% State vector:

% terrainState - terrain vertical position

% Observations:

% optRateX,  optRateY - line of sight (LOS) angular rate measurements rel 
% to sensor frame from a downwards looking optical flow sensor measured in
% rad/sec about the X and Y sensor axes. These rates are motion compensated.

% A positive LOS X rate is a RH rotation of the image about the X sensor
% axis, and is produced by either a positive ground relative velocity in 
% the direction of the Y axis.

% A positive LOS Y rate is a RH rotation of the image about the Y sensor
% axis, and is produced by either a negative ground relative velocity in 
% the direction of the X axis.

% range - laser range measurement from a sensor aligned with the Z body 
% axis and assuming a flat earth model

% Time varying parameters:

% q0,q1,q2,q3 - quaternion parameters defining the rotation from navigation 
% to body axes
% vel_x vel_y vel_z - NED flight vehicle velocities
% navPosD - vehicle vertical position

clear all;

%% define symbolic variables and constants
syms vel_x vel_y vel_z real % NED velocity : m/sec
syms R_LOS real % variance of LOS angular rate mesurements : (rad/sec)^2
syms R_RNG real % variance of range finder measurement : m^2
syms stateNoiseVar real % state process noise variance
syms navPosD real % position of vehicle in down axis : (m)
syms terrainState real % position of terrain in down axis : (m)
syms q0 q1 q2 q3 real % quaternions defining attitude of body axes relative to local NED
syms Popt real % state variance
nStates = 1;

%% define the process equations

% define the state update equation
terrainStateNew = terrainState;

%% derive the state transition matrix

% derive the state transition matrix
F = 1;

%% derive the covariance prediction equation

nextPopt = Popt + stateNoiseVar;

%% derive equations for sequential fusion of optical flow measurements

% derive the body to nav direction cosine matrix
Tbn = Quat2Tbn([q0,q1,q2,q3]);

% calculate relative velocity in sensor frame
relVelSensor = transpose(Tbn)*[vel_x;vel_y;vel_z];

% divide velocity by range to get predicted motion compensated flow rates
range = ((terrainState - navPosD)/Tbn(3,3));
optRateX =  relVelSensor(2)/range;
optRateY = -relVelSensor(1)/range;
predFlow = [optRateX;optRateY];
matlabFunction(predFlow,'file','calcPredFlow.m');

% calculate the observation jacobian
H_OPTX = jacobian(optRateX,terrainState);
H_OPTY = jacobian(optRateY,terrainState);
H_OPT = [H_OPTX;H_OPTY];
matlabFunction(H_OPT,'file','calcH_OPT.m');

% calculate Kalman gain vectors
% combine into a single K matrix to enable common expressions to be found
% note this matrix cannot be used in a single step fusion
K_OPTX = Popt*H_OPT(1)/(H_OPT(1)*Popt*H_OPT(1) + R_LOS);
K_OPTY = Popt*H_OPT(2)/(H_OPT(2)*Popt*H_OPT(2) + R_LOS);
K_OPT = [K_OPTX,K_OPTY];
matlabFunction(K_OPT,'file','calcK_OPT.m');

%% derive equations for fusion of range finder measurements
% assume laser is aligned in the X-Z body plane 
% calculate range from plane to centre of sensor fov assuming flat earth
range = ((terrainState - navPosD)/Tbn(3,3));
H_RNG = jacobian(range,terrainState); % measurement Jacobian
matlabFunction(H_RNG,'file','calcH_RNG.m');

% calculate the Kalman gain matrix and optimise algebra
K_RNG = (Popt*transpose(H_RNG))/(H_RNG*Popt*transpose(H_RNG) + R_RNG);
matlabFunction(K_RNG,'file','calcK_RNG.m');

%% derive the covariance update equations
%P = (I - K*H)*P
nextPX = Popt - K_OPT(1)*H_OPT(1)*Popt;
nextPY = Popt - K_OPT(2)*H_OPT(2)*Popt;
nextPR = Popt - K_RNG*H_RNG*Popt;

%% Save output
fileName = strcat('SymbolicOutput.mat');
save(fileName)