clear all;
load('2013-12-31 13-44.mat'); % first flight
%% IMU Data

zerotimestamp = IMU(1,1);
IMU(:,1) = (IMU(:,1) - zerotimestamp);
IMUtimestamp = IMU(:,1);
IMUtime  = IMU(:,2)*1e-3;
angRate  = IMU(:,3:5);
accel    = IMU(:,6:8);

%% GPS Data

if ~isempty(GPS)
    GPS(:,1) = (GPS(:,1) - zerotimestamp);
    goodDataIndices=find(GPS(:,2) == 3);
    GPStimestamp = GPS(goodDataIndices,1);
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
    MAG(:,1) = (MAG(:,1) - zerotimestamp);
    MAGtimestamp = MAG(:,1);
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
    NTUN(:,1) = (NTUN(:,1) - zerotimestamp);
    ADStimestamp = NTUN(:,1);
    ADStime  = NTUN(:,2)*1e-3;
    Veas     = NTUN(:,8);
    HgtBaro  = NTUN(:,9);
end

%% Reference Data

if ~isempty(ATT)
    ATT(:,1) = (ATT(:,1) - zerotimestamp);
    ATTtimestamp = ATT(:,1);
    ATTtime  = ATT(:,2)*1e-3;
    Roll     = ATT(:,3)*deg2rad;
    Pitch    = ATT(:,4)*deg2rad;
    Yaw      = ATT(:,5)*deg2rad;
end

%% Save to files
startTime = IMUtime(1,1);
alignTime = ATTtime(1)+1;
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

IMU = IMU;
GPS = GPS;
MAG = MAG;
ATT = ATT;
NTUN = NTUN;
save('IMU.txt','IMU','-ascii');
save('GPS.txt','GPS','-ascii');
save('MAG.txt','MAG','-ascii')
save('ATT.txt','ATT','-ascii');
save('NTUN.txt','NTUN','-ascii');
save('timing.txt','alignTime','startTime','endTime','msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay','EAS2TAS','-ascii');

clear all;
load('NavFilterTestData.mat');