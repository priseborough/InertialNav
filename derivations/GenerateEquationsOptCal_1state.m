% IMPORTANT - This script requires the Matlab symbolic toolbox
% Derivation of EKF equations for estimation of PX4Flow sensor errors
% Author:  Paul Riseborough
% Last Modified: 15 July 2014

% State vector:

% Optical flow sensor scale factor (allows for errors in focal length)

% Observations:

% line of sight (LOS) angular rate measurements (rel to sensor frame)
% from a downwards looking optical flow sensor measured in rad/sec about
% the X and Y sensor axes.

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
syms omegaX omegaY omegaZ real % bias corrected body rates - rad/sec
syms velN velE velD real % NED velocity - m/sec
syms R_LOS real % variance of LOS angular rate mesurements (rad/sec)^2
% syms a1 a2 a3 real % misalignment of the optical flow sensor (rad)
% syms w1_b w2_b real % bias errors on optical flow sensor rates
syms Kopt real % optical flow sensor LOS rate scale factor
syms range real % range estimate
%% define the process equations

% Define the state vector & number of states
%stateVector = [Kopt;a1;a2;a3;w1_b;w2_b];
stateVector = Kopt;
nStates=numel(stateVector);

% define the gyro bias state update equation
% w1_bNew = w1_b;
% w2_bNew = w2_b;

% define the misliagnment state update equation
% a1New = a1;
% a2New = a2;
% a3New = a3;

% define the optical flow scale factor update equation
KoptNew = Kopt;

% Define the process equations output vector
%processEqns = [KoptNew;a1New;a2New;a3New;w1_bNew;w2_bNew];
processEqns = KoptNew;

%% derive equations for sequential fusion of optical flow measurements

% define a symbolic covariance matrix using strings to represent
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

% define a symbolic nav to body direction cosine matrix using strings to represent
% '_l_' to represent '( '
% '_c_' to represent ,
% '_r_' to represent ')'
% these can be substituted later to create executable code
for rowIndex = 1:3
    for colIndex = 1:3
        eval(['syms Tnb_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['Tnb(',num2str(rowIndex),',',num2str(colIndex), ') = Tnb_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

% assume senosr axes are aligned with body axis plus a misalignment
% defined by 3 small angles about X, Y and Z axis
% Tsb is rotation from sensor to body
% Tbs = [ 1  ,  a3  , -a2 ; ...
%     -a3 ,  1   ,  a1 ; ...
%     a2 , -a1  ,  1 ];
% Tns = Tbs*Tnb;

% calculate relative velocity in sensor frame
% relVelSensor = Tns*[velN;velE;velD];
relVelSensor = Tnb*[velN;velE;velD];

% calculate angular rates in sensor frame
% omegaSensor = Tbs*[omegaX;omegaY;omegaZ];

% divide velocity by range, subtract body rates and apply scale factor to
% get predicted sensed angular optical rates relative to X and Y sensor axes
% optRateX =  Kopt*( relVelSensor(2)/range - omegaSensor(1));
% optRateY =  Kopt*(-relVelSensor(1)/range - omegaSensor(2));
optRateX =  Kopt*( relVelSensor(2)/range + omegaX);
optRateY =  Kopt*(-relVelSensor(1)/range + omegaY);

% calculate the observation jacobian
H_LOS = jacobian([optRateX;optRateY],stateVector);
[H_LOS,SH_LOS]=OptimiseAlgebra(H_LOS,'SH_LOS');

% calculate Kalman gain vectors
% combine into a single K matrix to enable common expressions to be found
% note this matrix cannot be used in a single step fusion
K_LOSX = (P*transpose(H_LOS(1,:)))/(H_LOS(1,:)*P*transpose(H_LOS(1,:)) + R_LOS);
K_LOSY = (P*transpose(H_LOS(2,:)))/(H_LOS(2,:)*P*transpose(H_LOS(2,:)) + R_LOS);
K_LOS = [K_LOSX,K_LOSY];
[K_LOS,SK_LOS]=OptimiseAlgebra(K_LOS,'SK_LOS');

%% Save output
fileName = strcat('SymbolicOutput',int2str(nStates),'.mat');
save(fileName);

%% Write equations to text and convert to m and c code fragments
SaveScriptCode(nStates);
ConvertToM(nStates);
ConvertToC(nStates);