clear all;
close all;
%load('NEO.mat')
load('search_pattern_log002.mat');

mintime = 10;
maxtime = 100;

i = 1;
prune_start_time = IMU.data(1,1);

while (i < size(IMU.data, 1) && (IMU.data(i,1) - prune_start_time) < maxtime)
    i = i+1;
end

IMU.data = IMU.data(1:i,:);

i = 1;

while (i < size(GPS.data, 1) && (GPS.data(i,1) - prune_start_time) < maxtime)
    i = i+1;
end

GPS.data = GPS.data(1:i,:);

i = 1;

while (i < size(AIRS.data, 1) && (AIRS.data(i,1) - prune_start_time) < maxtime)
    i = i+1;
end

AIRS.data = AIRS.data(1:i,:);

i = 1;

while (i < size(ATT.data, 1) && (ATT.data(i,1) - prune_start_time) < maxtime)
    i = i+1;
end

ATT.data = ATT.data(1:i,:);

%% Prune data according to time

%% IMU Data

zerotimestamp = IMU.data(1,1);
IMU.data(:,1) = (IMU.data(:,1) - zerotimestamp) * 1e3;
IMUtimestamp = IMU.data(:,1);
IMUtime  = IMU.data(:,1)*1e-3;
angRate  = IMU.data(:,5:7);
accel    = IMU.data(:,2:4);
    
MAGtimestamp = IMUtimestamp;
MAGtime  = IMUtime;
MagX     = IMU.data(:,8);
MagY     = IMU.data(:,9);
MagZ     = IMU.data(:,10);
MagBiasX = zeros(size(MagX,1), 1);
MagBiasY = zeros(size(MagY,1), 1);
MagBiasZ = zeros(size(MagZ,1), 1);

%% GPS Data

if ~isempty(GPS)
    GPS.data(:,1) = (GPS.data(:,1) - zerotimestamp) * 1e3;
    goodDataIndices=find(GPS.data(:,3) == 3);
    GPSstatus = GPS.data(goodDataIndices, 3);
    GPStimestamp = GPS.data(goodDataIndices,1);
    GPStime = GPS.data(goodDataIndices,1) * 1e-3;
    LatDeg = GPS.data(goodDataIndices,6);
    LngDeg = GPS.data(goodDataIndices,7);
    Hgt = GPS.data(goodDataIndices,8);
    deg2rad = pi/180;
    CourseDeg = GPS.data(goodDataIndices,12);
    VelN = GPS.data(goodDataIndices,9);
    VelE = GPS.data(goodDataIndices,10);
    VelD = GPS.data(goodDataIndices,11);
    % Get size right
    GndSpd = VelD;
    
    for i=1:size(VelN,1)
        GndSpd(i,1) = sqrt(VelN(i,1) * VelN(i,1) + VelE(i,1) * VelE(i,1));
    end
end

%% Magnetometer Data

% if ~isempty(MAG)
%     MAG.data(:,1) = (MAG.data(:,1) - zerotimestamp)*1000;
%     MAGtimestamp = MAG.data(:,1);
%     MAGtime  = MAG.data(:,2)*1e-3;
%     MagX     = MAG.data(:,3) - MAG.data(:,6);
%     MagY     = MAG.data(:,4) - MAG.data(:,7);
%     MagZ     = MAG.data(:,5) - MAG.data(:,8);
%     MagBiasX = - MAG.data(:,6);
%     MagBiasY = - MAG.data(:,7);
%     MagBiasZ = - MAG.data(:,8);
% end

%% Airspeed Data

if ~isempty(AIRS)
    AIRS.data(:,1) = (AIRS.data(:,1) - zerotimestamp) * 1e3;
    ADStimestamp = AIRS.data(:,1);
    ADStime  = AIRS.data(:,1) * 1e-3;
    Veas     = AIRS.data(:,3);
    % get size right
    HgtBaro  = zeros(size(AIRS.data(:,3), 1), 1);
    
    % Find closest baro reading
    if ~isempty(SENS)
        
        % AIRS is shorter than SENS (lower rate)
        % so we know i is a lower bound for j
        for i=1:size(AIRS.data,1)
            j = i;
            while ((SENS.data(j,1) - zerotimestamp) < ADStime(i, 1) && (j < size(SENS.data, 1)))
                j = j+1;
            end
            
            % we got a close measurement
            HgtBaro(i,1) = SENS.data(j,3);
        end
        
    end
end

%% Reference Data

if ~isempty(ATT.data)
    ATT.data(:,1) = (ATT.data(:,1) - zerotimestamp) * 1e3;
    ATTtimestamp = ATT.data(:,1);
    ATTtime  = ATT.data(:,1)*1e-3;
    Roll     = ATT.data(:,2) ./ pi .* 180;
    Pitch    = ATT.data(:,3) ./ pi .* 180;
    Yaw      = ATT.data(:,4) ./ pi .* 180;
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

% {'timestamp' 'TimeMS' 'GyrX' 'GyrY' 'GyrZ' 'AccX' 'AccY' 'AccZ'}
IMU = [IMUtimestamp IMUtimestamp angRate accel];
% {'timestamp' 'Status' 'TimeMS' 'Week' 'NSats' 'HDop' 'Lat' 'Lng' 'RelAlt' 'Alt' 'Spd' 'GCrs' 'VZ' 'T'}
gps_zeros = zeros(size(Hgt, 1), 1);
GPS = [GPStimestamp GPSstatus GPStimestamp gps_zeros gps_zeros gps_zeros LatDeg LngDeg Hgt Hgt GndSpd CourseDeg VelD gps_zeros];
GPSraw = [GPStimestamp LatDeg LngDeg Hgt VelN VelE VelD];
% {'timestamp' 'TimeMS' 'MagX' 'MagY' 'MagZ' 'OfsX' 'OfsY' 'OfsZ'}
MAG = [MAGtimestamp MAGtimestamp MagX MagY MagZ MagBiasX MagBiasY MagBiasZ];
% {'timestamp' 'TimeMS' 'Roll' 'Pitch' 'Yaw' 'ErrorRP' 'ErrorYaw'}
ATT = [ATTtimestamp ATTtimestamp Roll Pitch Yaw zeros(size(Yaw, 1), 1) zeros(size(Yaw, 1), 1)];
% {'timestamp' 'TimeMS' 'Yaw' 'WpDist' 'TargBrg' 'NavBrg' 'AltErr' 'Arspd' 'Alt' 'GSpdCM'}
NTUN = [ADStimestamp ADStimestamp zeros(size(Veas, 1), 1) zeros(size(Veas, 1), 1) zeros(size(Veas, 1), 1) zeros(size(Veas, 1), 1) zeros(size(Veas, 1), 1) Veas HgtBaro zeros(size(Veas, 1), 1)];
save('IMU.txt','IMU','-ascii');
save('GPS.txt','GPS','-ascii');
save('GPSraw.txt','GPSraw','-ascii');
save('MAG.txt','MAG','-ascii')
save('ATT.txt','ATT','-ascii');
save('NTUN.txt','NTUN','-ascii');
save('timing.txt','alignTime','startTime','endTime','msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay','EAS2TAS','-ascii');

clear all;
load('NavFilterTestData.mat');
