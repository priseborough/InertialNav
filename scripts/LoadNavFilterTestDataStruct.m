clear all;
%load('NEO.mat')
load('LEA.mat');

%% IMU Data

IMUframe = IMU.data(:,1);
IMUtime  = IMU.data(:,2)*1e-3;
angRate  = IMU.data(:,3:5);
accel    = IMU.data(:,6:8);

%% GPS Data

if ~isempty(GPS)
    goodDataIndices=find(GPS.data(:,2) == 3);
    GPSframe = GPS.data(goodDataIndices,1);
    GPStime = GPS.data(goodDataIndices,3);
    LatDeg = GPS.data(goodDataIndices,7);
    LngDeg = GPS.data(goodDataIndices,8);
    Hgt = GPS.data(goodDataIndices,9);
    deg2rad = pi/180;
    CourseDeg = GPS.data(goodDataIndices,12);
    GndSpd = GPS.data(goodDataIndices,11);
    VelD = GPS.data(goodDataIndices,13);
end

%% Magnetometer Data

if ~isempty(MAG)
    MAGframe = MAG.data(:,1);
    MAGtime  = MAG.data(:,2)*1e-3;
    MagX     = MAG.data(:,3) - MAG.data(:,6);
    MagY     = MAG.data(:,4) - MAG.data(:,7);
    MagZ     = MAG.data(:,5) - MAG.data(:,8);
    MagBiasX = - MAG.data(:,6);
    MagBiasY = - MAG.data(:,7);
    MagBiasZ = - MAG.data(:,8);
end

%% Air Data

if ~isempty(NTUN)
    ADSframe = NTUN.data(:,1);
    ADStime  = NTUN.data(:,2)*1e-3;
    Veas     = NTUN.data(:,8);
    HgtBaro  = NTUN.data(:,9);
end

%% Reference Data

if ~isempty(ATT.data)
    ATTframe = ATT.data(:,1);
    ATTtime  = ATT.data(:,2)*1e-3;
    Roll     = ATT.data(:,3)*deg2rad;
    Pitch    = ATT.data(:,4)*deg2rad;
    Yaw      = ATT.data(:,5)*deg2rad;
end

%% Save to files
save('NavFilterTestData.mat', ...
    'IMUframe','IMUtime','angRate','accel', ...
    'GPSframe','GPStime','LatDeg','LngDeg','Hgt','CourseDeg','GndSpd','VelD', ...
    'MAGframe','MAGtime','MagX','MagY','MagZ','MagBiasX','MagBiasY','MagBiasZ', ...
    'ATTframe','ATTtime','Roll','Pitch','Yaw', ...
    'ADSframe','ADStime','Veas','HgtBaro');
IMU = IMU.data;
GPS = GPS.data;
MAG = MAG.data;
ATT = ATT.data;
NTUN = NTUN.data;
save('../code/IMU.txt','IMU','-ascii');
save('../code/GPS.txt','GPS','-ascii');
save('../code/MAG.txt','MAG','-ascii')
save('../code/ATT.txt','ATT','-ascii');
save('../code/NTUN.txt','NTUN','-ascii');
clear all;
load('NavFilterTestData.mat');