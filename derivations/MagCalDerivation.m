%% define symbolic variables and constants
syms q0 q1 q2 q3 real % quaternions defining attitude of body axes relative to local NED
syms mx my mz real; % XYZ body fixed magnetic field bias - sensor units
syms mn me md real; % NED earth fixed magnetic field components - sensor units
syms kx ky kz real; % scale factor on magnetic field states from mot
syms mot real; % motor current input
syms R_MAG real  % variance for magnetic flux measurements - (sensor units)^2
syms Q_MAG real % magnetic field state process noise - (sensor units)^2

%% define the process equations

% Define the state vector & number of states
stateVector = [mn;me;md;mx;my;mz;kx;ky;kz];
nStates=numel(stateVector);

% define the process equations
mnNew = mn;
meNew = me;
mdNew = md;
mxNew = mx;
myNew = my;
mzNew = mz;
kxNew = kx;
kyNew = ky;
kzNew = kz;

% derive the body to nav direction cosine matrix
Tbn = Quat2Tbn([q0,q1,q2,q3]);

% Define the process equations output vector
processEqns = [mn;me;md;mx;my;mz;kx;ky;kz];

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

%% derive equations for fusion of magnetic field measurement
% define the measurements
magMeas = Tbn'*[mn;me;md] + mot*[kx;ky;kz] + [mx;my;mz];

H_MAG = jacobian(magMeas,stateVector); % measurement Jacobian
[H_MAG,SH_MAG]=OptimiseAlgebra(H_MAG,'SH_MAG');

K_MX = (P*transpose(H_MAG(1,:)))/(H_MAG(1,:)*P*transpose(H_MAG(1,:)) + R_MAG); % Kalman gain vector
[K_MX,SK_MX]=OptimiseAlgebra(K_MX,'SK_MX');
K_MY = (P*transpose(H_MAG(2,:)))/(H_MAG(2,:)*P*transpose(H_MAG(2,:)) + R_MAG); % Kalman gain vector
[K_MY,SK_MY]=OptimiseAlgebra(K_MY,'SK_MY');
K_MZ = (P*transpose(H_MAG(3,:)))/(H_MAG(3,:)*P*transpose(H_MAG(3,:)) + R_MAG); % Kalman gain vector
[K_MZ,SK_MZ]=OptimiseAlgebra(K_MZ,'SK_MZ');

%% Save output and convert to m and c code fragments
fileName = strcat('SymbolicOutput',int2str(nStates),'.mat');
save(fileName);
SaveScriptCode(nStates);
ConvertToM(nStates);
ConvertToC(nStates);
