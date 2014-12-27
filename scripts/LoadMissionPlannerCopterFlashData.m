clear all;
load('ekfOptFlowCopter.mat');
%maxtime = 600000;
%% IMU Data

zerotimestamp = IMU(1,2);
IMU(:,2) = (IMU(:,2) - zerotimestamp);
IMU(:,1) = IMU(:,2);
accZbias = mean(IMU(1:50,8)) + 9.8;
IMU(:,8) = IMU(:,8) - accZbias;
IMUtimestamp = IMU(:,2);
IMUtime  = IMU(:,2);
angRate  = IMU(:,3:5);
accel    = IMU(:,6:8);


%% GPS Data

if ~isempty(GPS)
    % trim first sample
    GPS = GPS(2:length(GPS),:);
    GPS(:,14) = (GPS(:,14) - zerotimestamp);
    goodDataIndices=find(GPS(:,2) == 3);
    GPSstatus = GPS(goodDataIndices,2);
    GPStimestamp = GPS(goodDataIndices,14);
    LatDeg = GPS(goodDataIndices,7);
    LngDeg = GPS(goodDataIndices,8);
    Hgt = GPS(goodDataIndices,9);
    deg2rad = pi/180;
    CourseDeg = GPS(goodDataIndices,12);
    GndSpd = GPS(goodDataIndices,11);
    VelD = GPS(goodDataIndices,13);
    gps_zeros = zeros(size(Hgt, 1), 1);
    GPS = [GPStimestamp GPSstatus GPStimestamp gps_zeros gps_zeros gps_zeros LatDeg LngDeg Hgt Hgt GndSpd CourseDeg VelD gps_zeros];
end

%% Magnetometer Data

if ~isempty(MAG)
    MAG(:,2) = (MAG(:,2) - zerotimestamp);
    MAGtimestamp = MAG(:,2);
    MAGtime  = MAG(:,2);
    MagX     = MAG(:,3);
    MagY     = MAG(:,4);
    MagZ     = MAG(:,5);
    MagBiasX = 0*MAG(:,6);
    MagBiasY = 0*MAG(:,7);
    MagBiasZ = 0*MAG(:,8);
    MAG = [MAGtimestamp MAGtimestamp MagX MagY MagZ MagBiasX MagBiasY MagBiasZ];
    
end

%% Air Data

if ~isempty(BARO)
    BARO(:,2) = (BARO(:,2) - zerotimestamp);
    ADStimestamp = BARO(:,2);
    ADStime  = BARO(:,2);
    Veas     = zeros(size(BARO(:,1)));
    HgtBaro  = BARO(:,3);
    NTUN = [ADStimestamp ADStimestamp zeros(size(Veas, 1), 1) zeros(size(Veas, 1), 1) zeros(size(Veas, 1), 1) zeros(size(Veas, 1), 1) zeros(size(Veas, 1), 1) Veas HgtBaro zeros(size(Veas, 1), 1)];
end

%% Range Finder Data

if ~isempty(EKF5)
    DIST(:,1) = (EKF5(:,2) - zerotimestamp);
    DIST(:,2) = EKF5(:,10);
    DIST(:,3) = 1;
end

%% PX4Flow Data

if (~isempty(OF) && ~isempty(EKF5))
    OF(:,2) = (OF(:,2) - zerotimestamp);
    EKF5(:,2) = (EKF5(:,2) - zerotimestamp);
    % timestamp, flowx, flowy, distance, quality, sensor id, angRateX,
    % angRateY
    % resample range measurements to be at same time coordinates as OF data 
    newIndex = 1;
    clear range rangeTime;
    rangeTime(1) = EKF5(1,2);
    range(1) = EKF5(1,10);
    for index = 1:length(EKF5(:,2))-1
        if ((EKF5(index + 1,2) - EKF5(index,2)) > 10)
            rangeTime(newIndex) = EKF5(index + 1,2);
            range(newIndex) = EKF5(index + 1,10);
            newIndex = newIndex +1;
        end
    end
    range = interp1(rangeTime,range,OF(:,2),'nearest',0);
    FLOW_OUT = [OF(:,2) -OF(:,4) -OF(:,5) range OF(:,3) 77*ones(size(OF(:,1))) OF(:,6) OF(:,7)];
end

%% Reference Data

% if ~isempty(AHR2)
%     AHR2(:,2) = (AHR2(:,2) - zerotimestamp);
%     ATTtimestamp = AHR2(:,2);
%     ATTtime  = AHR2(:,2);
%     Roll     = AHR2(:,3);
%     Pitch    = AHR2(:,4);
%     Yaw      = AHR2(:,5);
%     for i = 1:length(Yaw)
%         if Yaw(i) > 180
%             Yaw(i) = Yaw(i) - 360;
%         end
%     end
%     ATT = [ATTtimestamp ATTtimestamp Roll Pitch Yaw zeros(size(Yaw, 1), 1) zeros(size(Yaw, 1), 1)];
% end
if ~isempty(ATT)
    ATT(:,2) = (ATT(:,2) - zerotimestamp);
    ATTtimestamp = ATT(:,2);
    ATTtime  = ATT(:,2);
    Roll     = ATT(:,4);
    Pitch    = ATT(:,6);
    Yaw      = ATT(:,8);
    for i = 1:length(Yaw)
        if Yaw(i) > 180
            Yaw(i) = Yaw(i) - 360;
        end
    end
    ATT = [ATTtimestamp ATTtimestamp Roll Pitch Yaw zeros(size(Yaw, 1), 1) zeros(size(Yaw, 1), 1)];
end

%% Prune data
i = 1;
prune_start_time = IMU(1,1) + 1000

% Go right to the end of logfile minus 1000 ms
maxtime = IMU(end,1) - 1000

if (exist('IMU', 'var') == 1)
    while (i < size(IMU, 1) && (IMU(i,1) - prune_start_time) < maxtime)
        i = i+1;
    end
    
    IMU = IMU(1:i,:);
    
    i = 1;
end

if (exist('GPS', 'var') == 1)
    while (i < size(GPS, 1) && (GPS(i,1) - prune_start_time) < maxtime)
        i = i+1;
    end
    
    GPS = GPS(1:i,:);
    
    i = 1;
end

if (exist('MAG', 'var') == 1)
    while (i < size(MAG, 1) && (MAG(i,1) - prune_start_time) < maxtime)
        i = i+1;
    end
    
    MAG = MAG(1:i,:);
end

i = 1;

if (exist('FLOW', 'var') == 1)
    while (i < size(FLOW, 1) && (FLOW(i,1) - prune_start_time) < maxtime)
        i = i+1;
    end
    
    FLOW = FLOW(1:i,:);
end

i = 1;

if (exist('NTUN', 'var') == 1)
    while (i < size(NTUN, 1) && (NTUN(i,1) - prune_start_time) < maxtime)
        i = i+1;
    end
    
    NTUN = NTUN(1:i,:);
end

i = 1;

if (exist('ATT', 'var') == 1)
    while (i < size(ATT, 1) && (ATT(i,1) - prune_start_time) < maxtime)
        i = i+1;
    end
    
    ATT = ATT(1:i,:);
    
end

%% Save to files
startTime = IMUtime(1,1);
alignTime = (ATTtime(1)+1000)*0.001;
endTime = (max(IMUtime)-1000)*0.001;
msecVelDelay = 230;
msecPosDelay = 210;
msecHgtDelay = 350;
msecMagDelay = 30;
msecTasDelay = 210;
EAS2TAS = 1.0;

save('NavFilterTestData.mat', ...
    'IMUtimestamp','IMUtime','angRate','accel', ...
    'GPStimestamp','LatDeg','LngDeg','Hgt','CourseDeg','GndSpd','VelD', ...
    'MAGtimestamp','MAGtime','MagX','MagY','MagZ','MagBiasX','MagBiasY','MagBiasZ', ...
    'ATTtimestamp','ATTtime','Roll','Pitch','Yaw', ...
    'ADStimestamp','ADStime','Veas','HgtBaro', ...
    'alignTime','startTime','endTime', ...
    'msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay', ...
    'EAS2TAS', 'FLOW_OUT');

save('IMU.txt','IMU','-ascii');
save('GPS.txt','GPS','-ascii');
save('MAG.txt','MAG','-ascii')
save('ATT.txt','ATT','-ascii');
save('NTUN.txt','NTUN','-ascii');
save('timing.txt','alignTime','startTime','endTime','msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay','EAS2TAS','-ascii');
save('FLOW.txt','FLOW_OUT','-ascii');
save('DIST.txt','DIST','-ascii');

clear all;
load('NavFilterTestData.mat');