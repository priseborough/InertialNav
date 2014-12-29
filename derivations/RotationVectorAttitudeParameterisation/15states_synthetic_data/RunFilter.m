%% Set initial conditions
clear all;
dt = 1/100;
duration = 100;
indexLimit = round(duration/dt);
% define a yaw rate used to calculate the truth pose
yawRateTruth = 0.0;
yawAngleTruth = 0.0;
% define the truth magnetic field vector
magEarthTruth = [0.3;0.1;-0.5];
declination = atan2(magEarthTruth(2),magEarthTruth(1));
% define the sensor bias
magMeasBias = 0.05*[randn;randn;randn];
gyroRateBias = 5*(pi/180)*[randn;randn;randn];
% create output logging arrays
statesLog = zeros(15,indexLimit);
quatLog   = zeros(4,indexLimit);
velInnovLog = zeros(3,indexLimit);
magInnovLog = zeros(3,indexLimit);
velInnovVarLog = velInnovLog;
magInnovVarLog = magInnovLog;
angErrLog = zeros(1,indexLimit);
% Use a random initial orientation
quatTruth = [rand;randn;randn;randn];
quatLength = sqrt(quatTruth(1)^2 + quatTruth(2)^2 + quatTruth(3)^2 + quatTruth(4)^2);
quatTruth = quatTruth / quatLength;
TbnTruth = Quat2Tbn(quatTruth);
% initialise the filter to level
states = zeros(15,1);
quat = [1;0;0;0];
Tbn = Quat2Tbn(quat);
% define the initial state uncertainty
Sigma_velNED = 0.5; % 1 sigma uncertainty in horizontal velocity components
Sigma_dAngBias  = 5*pi/180*dt; % 1 Sigma uncertainty in delta angle bias
Sigma_angErr = 1; % 1 Sigma uncertainty in angular misalignment
Sigma_magNED = 0.1;
Sigma_magXYZ = 0.05;
covariance   = diag([Sigma_angErr*[1;1;10];Sigma_velNED*[1;1;1];Sigma_dAngBias*[1;1;1];Sigma_magNED*[1;1;1];Sigma_magXYZ*[1;1;1]].^2);
%% Main Loop
inhibitMag = 1;
headingAligned=0;
for index = 1:indexLimit
    % calculate the truth yaw angle
    yawAngleTruth = yawAngleTruth + dt*yawRateTruth;
    % synthesise IMU measurements
    angRate = 0.01*[randn;randn;randn] + [0;0;yawRateTruth] + gyroRateBias;
    accel = 0.01*[randn;randn;randn] + transpose(TbnTruth)*[0;0;-9.81];
    % predict states
    [quat, states, Tbn, delAng, delVel]  = PredictStates(quat,states,angRate,accel,dt);
    statesLog(:,index) = states;
    quatLog(:,index) = quat;
    % predict covariance matrix
    covariance  = PredictCovariance(delAng,delVel,quat,states,covariance,inhibitMag,dt);
    % synthesise velocity measurements
    measVel = [0;0;0];
    % fuse velocity measurements
    [quat,states,angErr,covariance,velInnov,velInnovVar] = FuseVelocity(quat,states,covariance,measVel);
    velInnovLog(:,index)    = velInnov;
    velInnovVarLog(:,index) = velInnovVar;
    angErrLog(1,index)      = angErr;
    % synthesise magnetometer measurements from truth data and add
    % measurement bias
    magMea = transpose(TbnTruth)*magEarthTruth + magMeasBias;
    % align the yaw angle and mag states
    if (headingAligned == 0 && index >= 5000 && angErr < 1E-4 )
        [quat,states] = AlignHeading(quat,states,magMea,declination);
        headingAligned=1;
        quat
        quatTruth
    end
    % fuse magnetometer measurements
    if (headingAligned == 1)
        [quat,states,covariance,magInnov,magInnovVar] = FuseMagnetometer(quat,states,covariance,magMea,inhibitMag,Tbn);
        magInnovLog(:,index) = magInnov;
        magInnovVarLog(:,index) = magInnovVar;
    end
end

