%% Set initial conditions
clear all;
load('fltTest.mat');
dt = 1/50;
startTime = 0.001*(IMU(1,2));
stopTime = 0.001*(IMU(length(IMU),2));
indexLimit = length(IMU);
magIndexlimit = length(MAG);
statesLog = zeros(10,indexLimit);
quatLog   = zeros(5,indexLimit);
velInnovLog = zeros(4,indexLimit);
angErrLog = zeros(2,indexLimit);
decInnovLog = zeros(2,magIndexlimit);
velInnovVarLog = velInnovLog;
decInnovVarLog = decInnovLog;
% initialise the filter to level
quat = [1;0;0;0];
states = zeros(9,1);
Tbn = Quat2Tbn(quat);
% Set the expected declination
measDec = 0.18;
% define the state covariances with the exception of the quaternion covariances
Sigma_velNED = 0.5; % 1 sigma uncertainty in horizontal velocity components
Sigma_dAngBias  = 5*pi/180*dt; % 1 Sigma uncertainty in delta angle bias
Sigma_angErr = 1; % 1 Sigma uncertainty in angular misalignment (rad)
covariance   = single(diag([Sigma_angErr*[1;1;10];Sigma_velNED*[1;1;1];Sigma_dAngBias*[1;1;1]].^2));
%% Main Loop
magIndex = 1;
time = 0;
angErr = 0;
headingAligned = 0;
for index = 1:indexLimit
    time=time+dt;
    % read IMU measurements and correct rates using estimated bias
    angRate = IMU(index,3:5)' - states(7:9)./dt;
    accel = IMU(index,6:8)';
    % predict states
    [quat, states, Tbn, delAng, delVel]  = PredictStates(quat,states,angRate,accel,dt);
    statesLog(1,index) = time;
    statesLog(2:10,index) = states;
    quatLog(1,index) = time;
    quatLog(2:5,index) = quat;
    % predict covariance matrix
    covariance  = PredictCovariance(delAng,delVel,quat,states,covariance,dt);
    % fuse velocity measurements - use synthetic measurements
    measVel = [0;0;0];
    [quat,states,angErr,covariance,velInnov,velInnovVar] = FuseVelocity(quat,states,covariance,measVel);
    velInnovLog(1,index) = time;
    velInnovLog(2:4,index) = velInnov;
    velInnovVarLog(1,index) = time;
    velInnovVarLog(2:4,index) = velInnovVar;
    angErrLog(1,index) = time;
    angErrLog(2,index) = angErr;
    % read magnetometer measurements
    while ((MAG(magIndex,1) < IMU(index,1)) && (magIndex < magIndexlimit))
        magIndex = magIndex + 1;
        magBody = 0.001*MAG(magIndex,3:5)';
        if (index > 250 && headingAligned==0 && angErr < 1e-4)
            quat = AlignHeading(quat,magBody,measDec);
            headingAligned = 1;
        end
        % fuse magnetometer measurements if new data available and when tilt has settled
        if (headingAligned == 1)
            [quat,states,covariance,decInnov,decInnovVar] = FuseMagnetometer(quat,states,covariance,magBody,measDec,Tbn);
            decInnovLog(1,magIndex) = time;
            decInnovLog(2,magIndex) = decInnov;
            decInnovVarLog(1,magIndex) = time;
            decInnovVarLog(2,magIndex) = decInnovVar;
        end
    end
end