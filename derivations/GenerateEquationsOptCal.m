% IMPORTANT - This script requires the Matlab symbolic toolbox
% Derivation of EKF equations for estimation of PX4Flow sensor errors
% Author:  Paul Riseborough
% Last Modified: 17 July 2014

% State vector:

% Optical flow sensor scale factor (allows for errors in focal length)

% Observations:

% line of sight (LOS) angular rate measurements (rel to sensor frame)
% from a downwards looking optical flow sensor measured in rad/sec about
% the X and Y sensor axes. These rates are not motion compensated.

% A positive LOS X rate is a RH rotation of the image about the X sensor
% axis, and is produced by either a negative sensor rotation rate about X
% or a positive ground relative velocity in the direction of the Y axis.

% A positive LOS Y rate is a RH rotation of the image about the Y sensor
% axis, and is produced by either a negative sensor rotation rate about Y
% or a negative ground relative velocity in the direction of the X axis.

% Time varying parameters:

% Tnb matrix defining the rotation from navigation to body axes
% Range from sensor to terrain
% NED flight vehicle velocities
% XYZ flight vehicle angular rates after corection for bias errors

clear all;

%% define symbolic variables and constants
syms omega_x omega_y omega_z real % bias corrected body rates - rad/sec
syms vel_x vel_y vel_z real % NED velocity - m/sec
syms R_OPT real % variance of LOS angular rate mesurements (rad/sec)^2
syms fScaleFactor real % focal length scale factor
syms range real % range estimate
syms fScaleFactorVar real % focal length scale factor variance
syms q0 q1 q2 q3 real % quaternions defining attitude of body axes relative to local NED
%% define the process equations

% Define the state vector & number of states
%stateVector = [Kopt;a1;a2;a3;w1_b;w2_b];
stateVector = fScaleFactor;
nStates=1;

% define the optical flow scale factor update equation
fScaleFactorNew = fScaleFactor;

% Define the process equations output vector
%processEqns = [KoptNew;a1New;a2New;a3New;w1_bNew;w2_bNew];
processEqns = fScaleFactorNew;

%% derive equations for sequential fusion of optical flow measurements

% derive the body to nav direction cosine matrix
Tbn = Quat2Tbn([q0,q1,q2,q3]);

% calculate relative velocity in sensor frame
relVelSensor = transpose(Tbn)*[vel_x;vel_y;vel_z];

% divide velocity by range, subtract body rates and apply scale factor to
% get predicted sensed angular optical rates relative to X and Y sensor axes
optRateX =  fScaleFactor*( relVelSensor(2)/range) - omega_x;
optRateY =  fScaleFactor*(-relVelSensor(1)/range) - omega_y;

% calculate the observation jacobian
H_OPTX = jacobian(optRateX,fScaleFactor);
H_OPTY = jacobian(optRateY,fScaleFactor);
H_OPT = [H_OPTX;H_OPTY];
[H_OPT,SH_OPT]=OptimiseAlgebra(H_OPT,'SH_OPT');

% calculate Kalman gain vectors
% combine into a single K matrix to enable common expressions to be found
% note this matrix cannot be used in a single step fusion
K_OPTX = (fScaleFactorVar*H_OPT(1))/(H_OPT(1)*fScaleFactorVar*H_OPT(1) + R_OPT);
K_OPTY = (fScaleFactorVar*H_OPT(2))/(H_OPT(2)*fScaleFactorVar*H_OPT(2) + R_OPT);
K_OPT = [K_OPTX,K_OPTY];
[K_OPT,SK_OPT]=OptimiseAlgebra(K_OPT,'SK_OPT');

%% Save output
fileName = strcat('SymbolicOutput',int2str(nStates),'.mat');
save(fileName);

%% Write equations to text and convert to m and c code fragments
SaveScriptCode(nStates);
ConvertToM(nStates);
ConvertToC(nStates);
