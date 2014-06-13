% IMPORTANT - This script requires the Matlab symbolic toolbox and takes ~3 hours to run

% Derivation of Navigation EKF using a local NED earth Tangent Frame and 
% XYZ body fixed frame
% Sequential fusion of velocity and position measurements
% Fusion of true airspeed
% Sequential fusion of magnetic flux measurements
% sequential fusion of angular LOS rate measurements from optical flow
% sensor assumed to be aligned witht eh Z body axis plus a small
% misalingment
% 23 state architecture.
% IMU data is assumed to arrive at a constant rate with a time step of dt
% IMU delta angle and velocity data are used as time varying parameters,
% not observations

% Author:  Paul Riseborough
% Last Modified: 23 Mar 2013

% State vector:
% quaternions (q0, q1, q2, q3)
% Velocity - m/sec (North, East, Down)
% Position - m (North, East, Down)
% Delta Angle bias - rad (X,Y,Z)
% Delta Velocity bias - m/s (Z)
% Wind Vector  - m/sec (North,East)
% Earth Magnetic Field Vector - milligauss (North, East, Down)
% Body Magnetic Field Vector - milligauss (X,Y,Z)
% Estimated position of terrain alogn D axis (m) 

% Observations:
% NED velocity - m/s
% NED position - m
% True airspeed - m/s
% XYZ magnetic flux - milligauss
% XY line of sight angular rate measurements from a downwards looking optical flow sensor
% range to terrain measurements

% Time varying parameters:
% XYZ delta angle measurements in body axes - rad
% XYZ delta velocity measurements in body axes - m/sec

clear all;

%% define symbolic variables and constants
syms dax day daz real % IMU delta angle measurements in body axes - rad
syms dvx dvy dvz real % IMU delta velocity measurements in body axes - m/sec
syms q0 q1 q2 q3 real % quaternions defining attitude of body axes relative to local NED
syms vn ve vd real % NED velocity - m/sec
syms pn pe pd real % NED position - m
syms dax_b day_b daz_b real % delta angle bias - rad
syms dvz_b real % delta velocity bias - m/sec
syms dt real % IMU time step - sec
syms gn ge gd real % NED gravity vector - m/sec^2
syms daxCov dayCov dazCov dvxCov dvyCov dvzCov real; % IMU delta angle and delta velocity measurement variances
syms vwn vwe real; % NE wind velocity - m/sec
syms magX magY magZ real; % XYZ body fixed magnetic field measurements - milligauss
syms magN magE magD real; % NED earth fixed magnetic field components - milligauss
syms R_VN R_VE R_VD real % variances for NED velocity measurements - (m/sec)^2
syms R_PN R_PE R_PD real % variances for NED position measurements - m^2
syms R_TAS real  % variance for true airspeed measurement - (m/sec)^2
syms R_MAG real  % variance for magnetic flux measurements - milligauss^2
syms R_LOS real % variance of LOS angular rate mesurements (rad/sec)^2
syms R_RNG real % variance of laser range finder observations
syms ptd real % location of terrain in D axis
syms a1 a2 a3 real % misalignment of the optical flow sensor (rad)
syms a4 real % pitch alignment of the laser range finder (rad)
%% define the process equations

% Define the state vector & number of states
stateVector = [q0;q1;q2;q3;vn;ve;vd;pn;pe;pd;dax_b;day_b;daz_b;dvz_b;vwn;vwe;magN;magE;magD;magX;magY;magZ;ptd];
nStates=numel(stateVector);

% define the measured Delta angle and delta velocity vectors
da = [dax; day; daz];
dv = [dvx; dvy; dvz];

% define the delta angle and delta velocity bias errors
da_b = [dax_b; day_b; daz_b];
dv_b = [0; 0; dvz_b];

% derive the body to nav direction cosine matrix
Tbn = Quat2Tbn([q0,q1,q2,q3]);

% define the bias corrected delta angles and velocities
dAngCor = da - da_b;
dVelCor = dv - dv_b;

% define the quaternion rotation vector
quat = [q0;q1;q2;q3];

% define the attitude update equations
delQuat = [1;
    0.5*dAngCor(1);
    0.5*dAngCor(2);
    0.5*dAngCor(3);
    ];
qNew = QuatMult(quat,delQuat);

% define the velocity update equations
vNew = [vn;ve;vd] + [gn;ge;gd]*dt + Tbn*dVelCor;

% define the position update equations
pNew = [pn;pe;pd] + [vn;ve;vd]*dt;

% define the IMU bias error update equations
dabNew = [dax_b; day_b; daz_b];
dvbNew = dvz_b;

% define the wind velocity update equations
vwnNew = vwn;
vweNew = vwe;

% define the earth magnetic field update equations
magNnew = magN;
magEnew = magE;
magDnew = magD;

% define the body magnetic field update equations
magXnew = magX;
magYnew = magY;
magZnew = magZ;

% define the terrain position update equation
ptdNew = ptd;

% Define the process equations output vector
processEqns = [qNew;vNew;pNew;dabNew;dvbNew;vwnNew;vweNew;magNnew;magEnew;magDnew;magXnew;magYnew;magZnew;ptdNew];

%% derive the state transition matrix

% derive the state transition matrix
F = jacobian(processEqns, stateVector);
[F,SF]=OptimiseAlgebra(F,'SF');

%% derive the covariance prediction equation
% This reduces the number of floating point operations by a factor of 6 or
% more compared to using the standard matrix operations in code

% Define the control (disturbance) vector. Error growth in the inertial
% solution is assumed to be driven by 'noise' in the delta angles and
% velocities, after bias effects have been removed. This is OK becasue we
% have sensor bias accounted for in the state equations.
distVector = [da;dv];

% derive the control(disturbance) influence matrix
G = jacobian(processEqns, distVector);
[G,SG]=OptimiseAlgebra(G,'SG');

% derive the state error matrix
imuNoise = diag([daxCov dayCov dazCov dvxCov dvyCov dvzCov]);
Q = G*imuNoise*transpose(G);
[Q,SQ]=OptimiseAlgebra(Q,'SQ');

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

% Derive the predicted covariance matrix using the standard equation
PP = F*P*transpose(F) + Q;

% Collect common expressions to optimise processing
[PP,SPP]=OptimiseAlgebra(PP,'SPP');

%% derive equations for sequential fusion of velocity and position measurements
H_VN= jacobian(vn,stateVector); % measurement Jacobian
K_VN = (P*transpose(H_VN))/(H_VN*P*transpose(H_VN) + R_VN);

H_VE= jacobian(ve,stateVector); % measurement Jacobian
K_VE = (P*transpose(H_VE))/(H_VE*P*transpose(H_VE) + R_VE);

H_VD= jacobian(vd,stateVector); % measurement Jacobian
K_VD = (P*transpose(H_VD))/(H_VD*P*transpose(H_VD) + R_VD);

H_PN= jacobian(pn,stateVector); % measurement Jacobian
K_PN = (P*transpose(H_PN))/(H_PN*P*transpose(H_PN) + R_PN);

H_PE= jacobian(pe,stateVector); % measurement Jacobian
K_PE = (P*transpose(H_PE))/(H_PE*P*transpose(H_PE) + R_PE);

H_PD= jacobian(pd,stateVector); % measurement Jacobian
K_PD = (P*transpose(H_PD))/(H_PD*P*transpose(H_PD) + R_PD);

% combine into a single H and K matrix (note these matrices cannot be used
% for a single step fusion, so each row|column must be used in a separate
% fusion step
H_VP  = [H_VN;H_VE;H_VD;H_PN;H_PE;H_PD];
clear    H_VN H_VE H_VD H_PN H_PE H_PD
K_VP = [K_VN,K_VE,K_VD,K_PN,K_PE,K_PD];
clear   K_VN K_VE K_VD K_PN K_PE K_PD
[K_VP,SK_VP]=OptimiseAlgebra(K_VP,'SK_VP');

%% derive equations for fusion of true airspeed measurements
VtasPred = sqrt((vn-vwn)^2 + (ve-vwe)^2 + vd^2); % predicted measurement
H_TAS = jacobian(VtasPred,stateVector); % measurement Jacobian
[H_TAS,SH_TAS]=OptimiseAlgebra(H_TAS,'SH_TAS'); % optimise processing
K_TAS = (P*transpose(H_TAS))/(H_TAS*P*transpose(H_TAS) + R_TAS);[K_TAS,SK_TAS]=OptimiseAlgebra(K_TAS,'SK_TAS'); % Kalman gain vector

%% derive equations for fusion of magnetic field measurement
magMeas = transpose(Tbn)*[magN;magE;magD] + [magX;magY;magZ]; % predicted measurement
H_MAG = jacobian(magMeas,stateVector); % measurement Jacobian
[H_MAG,SH_MAG]=OptimiseAlgebra(H_MAG,'SH_MAG');

K_MX = (P*transpose(H_MAG(1,:)))/(H_MAG(1,:)*P*transpose(H_MAG(1,:)) + R_MAG); % Kalman gain vector
[K_MX,SK_MX]=OptimiseAlgebra(K_MX,'SK_MX');
K_MY = (P*transpose(H_MAG(2,:)))/(H_MAG(2,:)*P*transpose(H_MAG(2,:)) + R_MAG); % Kalman gain vector
[K_MY,SK_MY]=OptimiseAlgebra(K_MY,'SK_MY');
K_MZ = (P*transpose(H_MAG(3,:)))/(H_MAG(3,:)*P*transpose(H_MAG(3,:)) + R_MAG); % Kalman gain vector
[K_MZ,SK_MZ]=OptimiseAlgebra(K_MZ,'SK_MZ');

%% derive equations for sequential fusion of optical flow measurements

% assume camera is aligned with Z body axis plus a misalignment
% defined by 3 small angles about X, Y and Z body axis
% Tsb is rotation from sensor to body
Tsb = [ 1  , -a3  ,  a2 ; ...
        a3 ,  1   , -a1 ; ...
       -a2 ,  a1  ,  1 ];
Tsn = Tbn*Tsb; 
% calculate range from plane to centre of sensor fov assuming flat earth
range = ((ptd - pd)/Tsn(3,3));
% calculate relative velocity in body frame
relVelBody = transpose(Tsn)*[vn;ve;vd];
% divide by range to get predicted angular LOS rates relative to X and Y
% axes
losRateX =  relVelBody(2)/range;
losRateY = -relVelBody(1)/range;

H_LOS = jacobian([losRateX;losRateY],stateVector); % measurement Jacobian
[H_LOS,SH_LOS]=OptimiseAlgebra(H_LOS,'SH_LOS');

% combine into a single K matrix to enable common expressions to be found
% note this matrix cannot be used in a single step fusion
K_LOSX = (P*transpose(H_LOS(1,:)))/(H_LOS(1,:)*P*transpose(H_LOS(1,:)) + R_LOS); % Kalman gain vector
K_LOSY = (P*transpose(H_LOS(2,:)))/(H_LOS(2,:)*P*transpose(H_LOS(2,:)) + R_LOS); % Kalman gain vector
K_LOS = [K_LOSX,K_LOSY];
[K_LOS,SK_LOS]=OptimiseAlgebra(K_LOS,'SK_LOS');

%% derive equations for fusion of laser range finder measurement

% assume laser is aligned in the X-Z body plane plus an alignment
% defined by rotation about the Y body axis. This allows the sensor to be
% aligned to accomodate different aircraft pitch trims. We would normally
% want the sensor to be aligned so that it was close to vertical during the
% phases of flight where height accuracy is most important
% Tsb is rotation from sensor to body
Tsb = [ 1       ,  0 , sin(a4) ; ...
        0       ,  1 , 0       ; ...
       -sin(a4) ,  0 , 1      ];
Tsn = Tbn*Tsb; 
% calculate range from plane to centre of sensor fov assuming flat earth
range = ((ptd - pd)/Tsn(3,3));
H_RNG = jacobian(range,stateVector); % measurement Jacobian
[H_RNG,SH_RNG]=OptimiseAlgebra(H_RNG,'SH_RNG');

% calculate the Kalman gain matrix and optimise algebra
K_RNG = (P*transpose(H_RNG))/(H_RNG*P*transpose(H_RNG) + R_RNG);
[K_RNG,SK_RNG]=OptimiseAlgebra(K_RNG,'SK_RNG');

%% Save output and convert to m and c code fragments
nStates = length(PP);
fileName = strcat('SymbolicOutput',int2str(nStates),'.mat');
save(fileName);
SaveScriptCode(nStates);
ConvertToM(nStates);
ConvertToC(nStates);