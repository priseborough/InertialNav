clear all;
%load('TestLog1.mat')
%load('TestLog2.mat')
load('TestLog3.mat')
%% IMU Data

for i = 1:length(IMU)-1
    dAngX(i) = 0.5*(IMU(i+1,3)+IMU(i,3))*(IMU(i+1,2)-IMU(i,2))*1e-3;
    dAngY(i) = 0.5*(IMU(i+1,4)+IMU(i,4))*(IMU(i+1,2)-IMU(i,2))*1e-3;
    dAngZ(i) = 0.5*(IMU(i+1,5)+IMU(i,5))*(IMU(i+1,2)-IMU(i,2))*1e-3;
    dVelX(i) = 0.5*(IMU(i+1,6)+IMU(i,6))*(IMU(i+1,2)-IMU(i,2))*1e-3;
    dVelY(i) = 0.5*(IMU(i+1,7)+IMU(i,7))*(IMU(i+1,2)-IMU(i,2))*1e-3;
    dVelZ(i) = 0.5*(IMU(i+1,8)+IMU(i,8))*(IMU(i+1,2)-IMU(i,2))*1e-3;
    IMUframe(i) = IMU(i+1,1);
    IMUtime(i) = (IMU(i+1,2)-IMU(2,2))*1e-3;
end

%% GPS Data

if ~isempty(GPS)
    goodDataIndices=find(GPS(:,2) == 3);
    GPSframe = GPS(goodDataIndices,1);
    LatDeg = GPS(goodDataIndices,7);
    LngDeg = GPS(goodDataIndices,8);
    Hgt = GPS(goodDataIndices,9);
    deg2rad = pi/180;
    CourseDeg = GPS(goodDataIndices,12);
    GndSpd = GPS(goodDataIndices,11);
    VelD = GPS(goodDataIndices,13);
end

%% Magnetometer Data

if ~isempty(MAG)
    MAGframe = MAG(:,1);
    MagX     = MAG(:,3) - MAG(:,6);
    MagY     = MAG(:,4) - MAG(:,7);
    MagZ     = MAG(:,5) - MAG(:,8);
    MagBiasX = - MAG(:,6);
    MagBiasY = - MAG(:,7);
    MagBiasZ = - MAG(:,8);
end

%% Air Data

if ~isempty(NTUN)
    ADSframe = NTUN(:,1);
    Veas     = NTUN(:,8);
    HgtBaro  = NTUN(:,9);
end

%% Reference Data

if ~isempty(ATT)
    ATTframe = ATT(:,1);
    Roll     = ATT(:,3)*deg2rad;
    Pitch    = ATT(:,4)*deg2rad;
    Yaw      = ATT(:,5)*deg2rad;
end

%% Save to file
save('NavFilterTestData.mat','dAngX','dAngY','dAngZ','dVelX','dVelY','dVelZ','IMUframe','IMUtime','GPSframe',...
    'LatDeg','LngDeg','Hgt','CourseDeg','GndSpd','VelD',...
    'MAGframe','MagX','MagY','MagZ','MagBiasX','MagBiasY','MagBiasZ',...
    'ATTframe','Roll','Pitch','Yaw', ...
    'ADSframe','Veas','HgtBaro');

clear all;
load('NavFilterTestData.mat');