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
syms omega_x omega_y omega_z real % bias corrected body rates : rad/sec
syms vel_x vel_y vel_z real % NED velocity : m/sec
syms R_OPT real % variance of LOS angular rate mesurements : (rad/sec)^2
syms R_RNG real % variance of range finder measurement : m^2
syms fScale real % focal length scale factor
syms stateNoiseVar1 stateNoiseVar2 real % state process noise variance
syms pd real % position of vehicle in down axis : (m)
syms ptd real % position of terrain in down axis : (m)
syms q0 q1 q2 q3 real % quaternions defining attitude of body axes relative to local NED
nStates = 2;
% define a symbolic covariance matrix using strings to represent the
% covariance matrix
% '_l_' to represent '( '
% '_c_' to represent ,
% '_r_' to represent ')'
% these can be substituted later to create executable code
for rowIndex = 1:nStates
    for colIndex = 1:nStates
        eval(['syms OP_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['P(',num2str(rowIndex),',',num2str(colIndex), ') = OP_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end
% define a symbolic body to NED direction cosine matrix using strings to represent
% '_l_' to represent '( '
% '_c_' to represent ,
% '_r_' to represent ')'
% these can be substituted later to create executable code
for rowIndex = 1:3
    for colIndex = 1:3
        eval(['syms Tbn_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['Tbn(',num2str(rowIndex),',',num2str(colIndex), ') = Tbn_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

%% define the process equations

% Define the state vector & number of states
stateVector = [fScale;ptd];

% define the state update equatiosn
fScaleNew = fScale;
ptdNew = ptd;

% Define the process equations output vector
%processEqns = [KoptNew;a1New;a2New;a3New;w1_bNew;w2_bNew];
processEqns = [fScaleNew;ptdNew];

%% derive the state transition matrix

% derive the state transition matrix
F = eye(2,2);

%% derive the covariance prediction equation
% define the state noise matrix
Q = [stateNoiseVar1 0;0 stateNoiseVar1];

% Derive the predicted covariance matrix using the standard equation
nextP = F*P*transpose(F) + Q;

%% derive equations for sequential fusion of optical flow measurements

% derive the body to nav direction cosine matrix
Tbn = Quat2Tbn([q0,q1,q2,q3]);

% calculate relative velocity in sensor frame
relVelSensor = transpose(Tbn)*[vel_x;vel_y;vel_z];

% divide velocity by range, subtract body rates and apply scale factor to
% get predicted sensed angular optical rates relative to X and Y sensor axes
range = ((ptd - pd)/Tbn(3,3));
optRateX =  fScale*( relVelSensor(2)/range) - omega_x;
optRateY =  fScale*(-relVelSensor(1)/range) - omega_y;

% calculate the observation jacobian
H_OPTX = jacobian(optRateX,stateVector);
H_OPTY = jacobian(optRateY,stateVector);
H_OPT = [H_OPTX;H_OPTY];
[H_OPT,SH_OPT]=OptimiseAlgebra(H_OPT,'SH_OPT');

% calculate Kalman gain vectors
% combine into a single K matrix to enable common expressions to be found
% note this matrix cannot be used in a single step fusion
K_OPTX = (P*transpose(H_OPT(1,:)))/(H_OPT(1,:)*P*transpose(H_OPT(1,:)) + R_OPT);
K_OPTY = (P*transpose(H_OPT(2,:)))/(H_OPT(2,:)*P*transpose(H_OPT(2,:)) + R_OPT);
K_OPT = [K_OPTX,K_OPTY];
[K_OPT,SK_OPT]=OptimiseAlgebra(K_OPT,'SK_OPT');

%% derive equations for fusion of range finder measurements
% assume laser is aligned in the X-Z body plane 
% calculate range from plane to centre of sensor fov assuming flat earth
range = ((ptd - pd)/Tbn(3,3));
H_RNG = jacobian(range,stateVector); % measurement Jacobian

% calculate the Kalman gain matrix and optimise algebra
K_RNG = (P*transpose(H_RNG))/(H_RNG*P*transpose(H_RNG) + R_RNG);
[K_RNG,SK_RNG]=OptimiseAlgebra(K_RNG,'SK_RNG');

%% derive the covariance update equations
%P = (I - K*H)*P
nextPX = (eye(2,2) - K_OPT(:,1)*H_OPT(1,:))*P;
nextPY = (eye(2,2) - K_OPT(:,2)*H_OPT(2,:))*P;
nextPR = (eye(2,2) - K_RNG*H_RNG)*P;

%% Save output
fileName = strcat('SymbolicOutput',int2str(nStates),'.mat');
save(fileName);

%% Write equations to text and convert to m and c code fragments
SaveScriptCode2State;
ConvertToM(nStates);
ConvertToC(nStates);
