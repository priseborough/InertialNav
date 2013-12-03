% Derivation of Navigation EKF using a local NED earth Tangent Frame and 
% XYZ body fixed frame
% Sequential fusion of velocity and position measurements
% Fusion of true airspeed
% Sequential fusion of magnetic flux measurements
% 24 state architecture.
% IMU data is assumed to arrive at a constant rate with a time step of dt
% IMU delta angle and velocity data are used as time varying parameters,
% not observations

% Author:  Paul Riseborough
% Last Modified: 3 Dec 2013

% State vector:
% quaternions (q0, q1, q2, q3)
% Velocity - m/sec (North, East, Down)
% Position - m (North, East, Down)
% Delta Angle bias - rad (X,Y,Z)
% Delta Velocity Bias - m/sec (X,Y,Z)
% Wind Vector  - m/sec (North,East)
% Earth Magnetic Field Vector - milligauss (North, East, Down)
% Body Magnetic Field Vector - milligauss (X,Y,Z)

% Observations:
% NED velocity - m/s
% NED position - m
% True airspeed - m/s
% XYZ magnetic flux - milligauss

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
syms dvx_b dvy_b dvz_b real % delta velocity bias - m/sec
syms dt real % IMU time step - sec
syms gn ge gd real % NED gravity vector - m/sec^2
syms omn ome omd real; % Earth rotation vector in local NED axes rad/sec - rad/sec
syms daxCov dayCov dazCov dvxCov dvyCov dvzCov real; % IMU delta angle and delta velocity measurement variances
syms vwn vwe real; % NE wind velocity - m/sec
syms magX magY magZ real; % XYZ body fixed magnetic field measurements - milligauss
syms magN magE magD real; % NED earth fixed magnetic field components - milligauss
syms R_VN R_VE R_VD real % variances for NED velocity measurements - (m/sec)^2
syms R_PN R_PE R_PD real % variances for NED position measurements - m^2
syms R_TAS real  % variance for true airspeed measurement - (m/sec)^2
syms R_MAG real  % variance for magnetic flux measurements - milligauss^2

%% define the process equations

% Define the state vector & number of states
stateVector = [q0;q1;q2;q3;vn;ve;vd;pn;pe;pd;dax_b;day_b;daz_b;dvx_b;dvy_b;dvz_b;vwn;vwe;magN;magE;magD;magX;magY;magZ];
nStates=numel(stateVector);

% define the measured Delta angle and delta velocity vectors
da = [dax; day; daz];
dv = [dvx; dvy; dvz];

% define the delta angle and delta velocity bias errors
da_b = [dax_b; day_b; daz_b];
dv_b = [dvx_b; dvy_b; dvz_b];

% derive the body to nav direction cosine matrix
Tbn = Quat2Tbn([q0,q1,q2,q3]);

% define the bias corrected delta angle
% Ignore coning compensation and earths rotation as these effect are
% negligible in terms of covariance growth compared to other efects for our
% grade of sensor
% deltaAngle = da - da_b + 1/12*cross(da_prev,da) - transpose(Cbn)*([omn; ome; omd])*dt;
deltaAngle = da - da_b;

% define the bias corrected delta velocity
% Ignore sculling as this effect is negligible in terms of covariance growth 
% compared to other effects for our grade of sensor
% deltaVelocity = dv - dv_b + 0.5*cross(da,dv) + 1/12*(cross(da_prev,dv) + cross(dv_prev,da));
deltaVelocity = dv - dv_b;

% define the quaternion rotation vector
quat = [q0;q1;q2;q3];

% define the attitude update equations
% use a first order expansion of rotation to calculate the quaternion increment
% acceptable for propagation of covariances
delQuat = [1;
    0.5*deltaAngle(1);
    0.5*deltaAngle(2);
    0.5*deltaAngle(3);
    ];
qNew = QuatMult(quat,delQuat);

% define the velocity update equations
vNew = [vn;ve;vd] + [gn;ge;gd]*dt + Tbn*deltaVelocity;% - cross(2*[omn; ome; omd],[vn;ve;vd])*dt;

% define the position update equations
pNew = [pn;pe;pd] + [vn;ve;vd]*dt;

% define the IMU bias error update equations
dabNew = [dax_b; day_b; daz_b];
dvbNew = [dvx_b; dvy_b; dvz_b];

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

% Define the process equations output vector
processEqns = [qNew;vNew;pNew;dabNew;dvbNew;vwnNew;vweNew;magNnew;magEnew;magDnew;magXnew;magYnew;magZnew];

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
% '_lp_' to represent '( '
% '_c_' to represent ,
% '_rp_' to represent ')' 
% these can be substituted later to create executable code
for rowIndex = 1:nStates
    for colIndex = 1:nStates
        eval(['syms P_lp_',num2str(rowIndex),'_c_',num2str(colIndex), '_rp_ real']);
        eval(['P(',num2str(rowIndex),',',num2str(colIndex), ') = P_lp_',num2str(rowIndex),'_c_',num2str(colIndex),'_rp_;']);
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
% for a single step fusion, so each row|column mst be used in a separate
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

% combine into a single K matrix (note this matrix cannot be used
% for a single step fusion, so each column mst be used sequentially
K_MAG = [K_MX,K_MY,K_MZ];
clear    K_MX K_MY K_MZ

%% Open output file
fid = fopen('SymbolicOutput.m','wt');

%% Write equation for state transition matrix
fprintf(fid,'SF = zeros(%d,1);\n',numel(SF));
for rowIndex = 1:numel(SF)
    string = char(SF(rowIndex,1));
    fprintf(fid,'SF(%d) = %s;\n',rowIndex,string);
end

fprintf(fid,'\n');
fprintf(fid,'F = zeros(%d,%d);\n',nStates,nStates);
for rowIndex = 1:nStates
    for colIndex = 1:nStates
        string = char(F(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'F(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end

%% Write equations for control influence (disturbance) matrix
fprintf(fid,'\n');
fprintf(fid,'SG = zeros(%d,1);\n',numel(SG));
for rowIndex = 1:numel(SG)
    string = char(SG(rowIndex,1));
    fprintf(fid,'SG(%d) = %s;\n',rowIndex,string);
end

fprintf(fid,'\n');
fprintf(fid,'G = zeros(%d,%d);\n',nStates,numel([da;dv]));
for rowIndex = 1:nStates
    for colIndex = 1:numel([da;dv])
        string = char(G(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'G(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end

%% Write equations for state error matrix
fprintf(fid,'\n');
fprintf(fid,'SQ = zeros(%d,1);\n',numel(SQ));
for rowIndex = 1:numel(SQ)
    string = char(SQ(rowIndex,1));
    fprintf(fid,'SQ(%d) = %s;\n',rowIndex,string);
end

fprintf(fid,'\n');
fprintf(fid,'Q = zeros(%d,%d);\n',nStates,nStates);
for rowIndex = 1:nStates
    for colIndex = 1:nStates
        string = char(Q(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'Q(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end

%% Write equations for covariance prediction
fprintf(fid,'\n');
fprintf(fid,'SPP = zeros(%d,1);\n',numel(SPP));
for rowIndex = 1:numel(SPP)
    string = char(SPP(rowIndex,1));
    fprintf(fid,'SPP(%d) = %s;\n',rowIndex,string);
end

fprintf(fid,'\n');
fprintf(fid,'PP = zeros(%d,%d);\n',nStates,nStates);
for rowIndex = 1:nStates
    for colIndex = 1:nStates
        string = char(PP(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'PP(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end

%% Write equations for velocity and position data fusion
[nRow,nCol] = size(H_VP);
fprintf(fid,'\n');
fprintf(fid,'H_VP = zeros(%d,%d);\n',nRow,nCol);
for rowIndex = 1:nRow
    for colIndex = 1:nCol
        string = char(H_VP(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'H_VP(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end

[nRow,nCol] = size(SK_VP);
fprintf(fid,'\n');
fprintf(fid,'SK_VP = zeros(%d,%d);\n',nRow,nCol);
for rowIndex = 1:nRow
    for colIndex = 1:nCol
        string = char(SK_VP(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'SK_VP(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end

[nRow,nCol] = size(K_VP);
fprintf(fid,'\n');
fprintf(fid,'K_VP = zeros(%d,%d);\n',nRow,nCol);
for rowIndex = 1:nRow
    for colIndex = 1:nCol
        string = char(K_VP(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'K_VP(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end

%% Write equations for true airspeed data fusion
fprintf(fid,'\n');
fprintf(fid,'SH_TAS = zeros(%d,1);\n',numel(SH_TAS));
for rowIndex = 1:numel(SH_TAS)
    string = char(SH_TAS(rowIndex,1));
    fprintf(fid,'SH_TAS(%d) = %s;\n',rowIndex,string);
end

[nRow,nCol] = size(H_TAS);
fprintf(fid,'\n');
fprintf(fid,'H_TAS = zeros(%d,%d);\n',nRow,nCol);
for rowIndex = 1:nRow
    for colIndex = 1:nCol
        string = char(H_TAS(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'H_TAS(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end

fprintf(fid,'\n');
fprintf(fid,'SK_TAS = zeros(%d,1);\n',numel(SK_TAS));
for rowIndex = 1:numel(SK_TAS)
    string = char(SK_TAS(rowIndex,1));
    fprintf(fid,'SK_TAS(%d) = %s;\n',rowIndex,string);
end

[nRow,nCol] = size(K_TAS);
fprintf(fid,'\n');
fprintf(fid,'K_TAS = zeros(%d,%d);\n',nRow,nCol);
for rowIndex = 1:nRow
    for colIndex = 1:nCol
        string = char(K_TAS(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'K_TAS(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end

%% Write equations for magnetometer data fusion
fprintf(fid,'\n');
fprintf(fid,'SH_MAG = zeros(%d,1);\n',numel(SH_MAG));
for rowIndex = 1:numel(SH_MAG)
    string = char(SH_MAG(rowIndex,1));
    fprintf(fid,'SH_MAG(%d) = %s;\n',rowIndex,string);
end

[nRow,nCol] = size(H_MAG);
fprintf(fid,'\n');
fprintf(fid,'H_MAG = zeros(%d,%d);\n',nRow,nCol);
for rowIndex = 1:nRow
    for colIndex = 1:nCol
        string = char(H_MAG(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'H_MAG(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end

fprintf(fid,'\n');
fprintf(fid,'SK_MX = zeros(%d,1);\n',numel(SK_MX));
for rowIndex = 1:numel(SK_MX)
    string = char(SK_MX(rowIndex,1));
    fprintf(fid,'SK_MX(%d) = %s;\n',rowIndex,string);
end

fprintf(fid,'\n');
fprintf(fid,'SK_MY = zeros(%d,1);\n',numel(SK_MY));
for rowIndex = 1:numel(SK_MY)
    string = char(SK_MY(rowIndex,1));
    fprintf(fid,'SK_MY(%d) = %s;\n',rowIndex,string);
end

fprintf(fid,'\n');
fprintf(fid,'SK_MZ = zeros(%d,1);\n',numel(SK_MZ));
for rowIndex = 1:numel(SK_MZ)
    string = char(SK_MZ(rowIndex,1));
    fprintf(fid,'SK_MZ(%d) = %s;\n',rowIndex,string);
end

[nRow,nCol] = size(K_MAG);
fprintf(fid,'\n');
fprintf(fid,'K_MAG = zeros(%d,%d);\n',nRow,nCol);
for rowIndex = 1:nRow
    for colIndex = 1:nCol
        string = char(K_MAG(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'K_MAG(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end

%% Close output file
fclose(fid);