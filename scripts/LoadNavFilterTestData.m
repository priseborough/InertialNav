clear all;
%load('TestLog1.mat')
%load('TestLog2.mat')
%load('TestLog3.mat')

%% IMU Data

IMUframe = IMU(:,1);
IMUtime  = IMU(:,2)*1e-3;
angRate  = IMU(:,3:5);
accel    = IMU(:,6:8);

%% GPS Data

if ~isempty(GPS)
    goodDataIndices=find(GPS(:,2) == 3);
    GPSframe = GPS(goodDataIndices,1);
    GPStime = GPS(goodDataIndices,3);
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
    MAGtime  = MAG(:,2)*1e-3;
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
    ADStime  = NTUN(:,2)*1e-3;
    Veas     = NTUN(:,8);
    HgtBaro  = NTUN(:,9);
end

%% Reference Data

if ~isempty(ATT)
    ATTframe = ATT(:,1);
    ATTtime  = ATT(:,2)*1e-3;
    Roll     = ATT(:,3)*deg2rad;
    Pitch    = ATT(:,4)*deg2rad;
    Yaw      = ATT(:,5)*deg2rad;
end

%% Save to files
save('NavFilterTestData.mat', ...
    'IMUframe','IMUtime','angRate','accel', ...
    'GPSframe','GPStime','LatDeg','LngDeg','Hgt','CourseDeg','GndSpd','VelD', ...
    'MAGframe','MAGtime','MagX','MagY','MagZ','MagBiasX','MagBiasY','MagBiasZ', ...
    'ATTframe','ATTtime','Roll','Pitch','Yaw', ...
    'ADSframe','ADStime','Veas','HgtBaro');
save('../code/IMU.txt','IMU','-ascii');
save('../code/GPS.txt','GPS','-ascii');
save('../code/MAG.txt','MAG','-ascii')
save('../code/ATT.txt','ATT','-ascii');
save('../code/NTUN.txt','NTUN','-ascii');
clear all;
load('NavFilterTestData.mat');