clear all;
%load('NEO.mat')
%load('LEA.mat');
%load('SITL1.mat'); % sim_gps_delay = 0
%load('SITL2.mat'); % sim_gps_delay = 3
%load('SITL3.mat'); % sim_gps_delay = 2
%load('SITL4.mat'); % sim_gps_delay = 1
load('test2.mat');

%% IMU Data

zerotimestamp = IMU.data(1,1);
IMU.data(:,1) = (IMU.data(:,1) - zerotimestamp)*1000;
IMUtimestamp = IMU.data(:,1);
IMUtime  = IMU.data(:,2)*1e-3;
angRate  = IMU.data(:,3:5);
accel    = IMU.data(:,6:8);

%% GPS Data

if ~isempty(GPS)
    GPS.data(:,1) = (GPS.data(:,1) - zerotimestamp)*1000;
    goodDataIndices=find(GPS.data(:,2) == 3);
    GPStimestamp = GPS.data(goodDataIndices,1);
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
    MAG.data(:,1) = (MAG.data(:,1) - zerotimestamp)*1000;
    MAGtimestamp = MAG.data(:,1);
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
    NTUN.data(:,1) = (NTUN.data(:,1) - zerotimestamp)*1000;
    ADStimestamp = NTUN.data(:,1);
    ADStime  = NTUN.data(:,2)*1e-3;
    Veas     = NTUN.data(:,8);
    HgtBaro  = NTUN.data(:,9);
end

%% Reference Data

if ~isempty(ATT.data)
    ATT.data(:,1) = (ATT.data(:,1) - zerotimestamp)*1000;
    ATTtimestamp = ATT.data(:,1);
    ATTtime  = ATT.data(:,2)*1e-3;
    Roll     = ATT.data(:,3)*deg2rad;
    Pitch    = ATT.data(:,4)*deg2rad;
    Yaw      = ATT.data(:,5)*deg2rad;
end

%% Save to files
% alignTime = min(IMUtime(IMUtimestamp>GPStimestamp(find(GndSpd  >8, 1 )))) - 10;
% startTime = max((alignTime - 30),(IMUtime(1,1) + 1));
% alignTime = max((startTime+1),alignTime);
startTime = (IMUtime(1,1) + 1);
alignTime = startTime + 10;
alignTime = max((startTime+1),alignTime);
endTime = max(IMUtime)-10;
msecVelDelay = 230;
msecPosDelay = 210;
msecHgtDelay = 350;
msecMagDelay = 30;
msecTasDelay = 210;
EAS2TAS = 1.0;

save('NavFilterTestData.mat', ...
    'IMUtimestamp','IMUtime','angRate','accel', ...
    'GPStimestamp','GPStime','LatDeg','LngDeg','Hgt','CourseDeg','GndSpd','VelD', ...
    'MAGtimestamp','MAGtime','MagX','MagY','MagZ','MagBiasX','MagBiasY','MagBiasZ', ...
    'ATTtimestamp','ATTtime','Roll','Pitch','Yaw', ...
    'ADStimestamp','ADStime','Veas','HgtBaro', ...
    'alignTime','startTime','endTime', ...
    'msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay', ...
    'EAS2TAS');

IMU = IMU.data;
GPS = GPS.data;
MAG = MAG.data;
ATT = ATT.data;
NTUN = NTUN.data;
save('IMU.txt','IMU','-ascii');
save('GPS.txt','GPS','-ascii');
save('MAG.txt','MAG','-ascii')
save('ATT.txt','ATT','-ascii');
save('NTUN.txt','NTUN','-ascii');
save('timing.txt','alignTime','startTime','endTime','msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay','EAS2TAS','-ascii');

clear all;
load('NavFilterTestData.mat');