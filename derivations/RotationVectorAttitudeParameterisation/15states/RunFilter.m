%% Set initial conditions
clear all;
dt = 1/100;
duration = 100;
indexLimit = round(duration/dt);
% define a yaw rate used to calculate the truth pose
yawRateTruth = 0.5;
yawAngleTruth = 0.0;
% define the truth magnetic field vector
magEarthTruth = [0.3;0.1;-0.5];
% create output logging arrays
statesLog = zeros(15,indexLimit);
quatLog   = zeros(4,indexLimit);
velInnovLog = zeros(3,indexLimit);
magInnovLog = zeros(3,indexLimit);
velInnovVarLog = velInnovLog;
magInnovVarLog = magInnovLog;
% Initialise the orientation and earth field estimate to bad values
states = zeros(15,1);
%quat = [1/sqrt(2);1/sqrt(2);0;0];
quat = [rand;randn;randn;randn];
quatMag = sqrt(quat(1)^2 + quat(2)^2 + quat(3)^2 + quat(4)^2);
quat(1:4) = quat / quatMag;
Tbn = Quat2Tbn(quat);
magMeasBias = 0.05*[randn;randn;randn];
% initialise the earth field using the initial orientation and measurements
states(10:12) = Tbn*(transpose(Tbn)*(magEarthTruth) + magMeasBias);
% define the state covariances with the exception of the quaternion covariances
Sigma_velNED = 0.5; % 1 sigma uncertainty in horizontal velocity components
Sigma_dAngBias  = 5*pi/180*dt; % 1 Sigma uncertainty in delta angle bias
Sigma_angErr = 10*pi/180; % 1 Sigma uncertainty in angular misalignment
Sigma_magNED = 0.5;
Sigma_magXYZ = 0.05;
covariance   = diag([Sigma_angErr*[1;1;1];Sigma_velNED*[1;1;1];Sigma_dAngBias*[1;1;1];Sigma_magNED*[1;1;1];Sigma_magXYZ*[1;1;1]].^2);
% define the IMU bias errors
gyroRateBias = 5*(pi/180)*[randn;randn;randn];
%% Main Loop
for index = 1:indexLimit
    % calculate the truth yaw angle
    yawAngleTruth = yawAngleTruth + dt*yawRateTruth;
    % synthesise IMU measurements
    angRate = 0.01*[randn;randn;randn] + [0;0;yawRateTruth] + gyroRateBias;
    accel = 0.01*[randn;randn;randn] + [0;0;-9.81];
    % predict states
    [quat, states, Tbn, delAng, delVel]  = PredictStates(quat,states,angRate,accel,dt);
    statesLog(:,index) = states;
    quatLog(:,index) = quat;
    % predict covariance matrix
    covariance  = PredictCovariance(delAng,delVel,quat,states,covariance,dt);
    % synthesise magnetometer measurements from truth data and add
    % measurement bias
    magMea = [magEarthTruth(1)*cos(yawAngleTruth)+magEarthTruth(2)*sin(yawAngleTruth);magEarthTruth(2)*cos(yawAngleTruth)-magEarthTruth(1)*sin(yawAngleTruth);magEarthTruth(3)];
    magMea = magMea + magMeasBias;
    % fuse magnetometer measurements
    [quat,states,covariance,magInnov,magInnovVar] = FuseMagnetometer(quat,states,covariance,magMea,Tbn);
    magInnovLog(:,index) = magInnov;
    magInnovVarLog(:,index) = magInnovVar;
    % synthesise velocity measurements
    measVel = [0;0;0];
    % fuse velocity measurements
    [quat,states,covariance,velInnov,velInnovVar] = FuseVelocity(quat,states,covariance,measVel);
    velInnovLog(:,index) = velInnov;
    velInnovVarLog(:,index) = velInnovVar;
end

