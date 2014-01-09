dirName = '../plots/';
t0 = EKF1.data(1,1);
%% Euler Angles
fileName = 'EulerAngleEstimates';
figure;
subplot(3,1,1);
if (exist('AHR2','var'))
    plot(AHR2.data(:,1)-t0,AHR2.data(:,3),'g');
    hold on;
end
plot(EKF1.data(:,1)-t0,0.01*EKF1.data(:,3),'b');
hold off;
grid on;
ylim([-180 180]);
xlabel('time (sec)');ylabel('roll (deg)');
title('Euler Angle Estimates');
subplot(3,1,2);
if (exist('AHR2','var'))
    plot(AHR2.data(:,1)-t0,AHR2.data(:,4),'g');
    hold on;
end
plot(EKF1.data(:,1)-t0,0.01*EKF1.data(:,4),'b');
hold off;
grid on;
ylim([-180 180]);
xlabel('time (sec)');ylabel('pitch (deg)');
subplot(3,1,3);
if (exist('AHR2','var'))
    plot(AHR2.data(:,1)-t0,AHR2.data(:,5),'g');
    hold on;
end
plot(EKF1.data(:,1)-t0,0.01*EKF1.data(:,5),'b');
hold off;
grid on;
ylim([0 360]);
xlabel('time (sec)');ylabel('yaw (deg)');
saveas(gcf,strcat(dirName,fileName,'.jpg'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% NED velocity
fileName = 'VelocityEstimates';
figure;
subplot(3,1,1);
if (exist('GPS','var'))
    plot(GPS.data(:,1)-t0,GPS.data(:,11).*cos(GPS.data(:,12)*pi/180),'g');
    hold on;
end
plot(EKF1.data(:,1)-t0,EKF1.data(:,6),'b');
hold off;
grid on;
xlabel('time (sec)');ylabel('North Velocity (m/s)');
title('NED Velocity Estimates');
subplot(3,1,2);
if (exist('GPS','var'))
    plot(GPS.data(:,1)-t0,GPS.data(:,11).*sin(GPS.data(:,12)*pi/180),'g');
    hold on;
end
plot(EKF1.data(:,1)-t0,EKF1.data(:,7),'b');
hold off;
hold off;
grid on;
xlabel('time (sec)');ylabel('East Velocity (m/s)');
subplot(3,1,3);
if (exist('GPS','var'))
    plot(GPS.data(:,1)-t0,GPS.data(:,13),'g');
    hold on;
end
plot(EKF1.data(:,1)-t0,EKF1.data(:,8),'b');
hold off;
grid on;
xlabel('time (sec)');ylabel('Down Velocity (m/s)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% NE Position and Height
if (exist('GPS','var'))
    earthRadius = 6378100;
    lat = GPS.data(:,7)*pi/180;
    lon = GPS.data(:,8)*pi/180;
    latRef = median(lat(1:10));
    lonRef = median(lon(1:10));
    posNgps = earthRadius * (lat - latRef);
    posEgps = earthRadius * cos(latRef) * (lon - lonRef);
    clear lat lon latRef lonRef;
end
fileName = 'PositionEstimates';
figure;
subplot(3,1,1);
if (exist('GPS','var'))
    plot(GPS.data(:,1)-t0,posNgps,'g');hold on;
end
plot(EKF1.data(:,1)-t0,EKF1.data(:,9),'b');hold off
grid on;
xlabel('time (sec)');ylabel('North Position (m)');
title('NED Position Estimates');
subplot(3,1,2);
if (exist('GPS','var'))
    plot(GPS.data(:,1)-t0,posEgps,'g');hold on;
end
plot(EKF1.data(:,1)-t0,EKF1.data(:,10),'b');hold off;
grid on;
xlabel('time (sec)');ylabel('East Position (m)');
subplot(3,1,3);
if (exist('NTUN','var'))
    plot(NTUN.data(:,1)-t0,NTUN.data(:,9)-NTUN.data(find(NTUN.data(:,9)~=0, 1 ),9),'g');hold on;
end
plot(EKF1.data(:,1)-t0,-EKF1.data(:,11),'b');hold off
grid on;
ylim([-10,max(-EKF1.data(:,11))+10]);
xlabel('time (sec)');ylabel('Height (m)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% Gyro Bias
fileName = 'GyroBiasEstimates';
figure;
subplot(3,1,1);
plot(EKF1.data(:,1)-t0,0.01*EKF1.data(:,12),'b');
grid on;
xlabel('time (sec)');ylabel('X bias (deg/min)');
title('Gyro Bias Error Estimates');
subplot(3,1,2);
plot(EKF1.data(:,1)-t0,0.01*EKF1.data(:,13),'b');
grid on;
xlabel('time (sec)');ylabel('Y bias (deg/min)');
subplot(3,1,3);
plot(EKF1.data(:,1)-t0,0.01*EKF1.data(:,14),'b');
grid on;
xlabel('time (sec)');ylabel('Z bias (deg/min)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% Accelerometer Bias
fileName = 'AccelBiasEstimates';
figure;
subplot(3,1,1);
plot(EKF2.data(:,1)-t0,EKF2.data(:,3),'b');
grid on;
xlabel('time (sec)');ylabel('X bias (cm/s^2)');
title('Accelerometer Bias Error Estimates');
subplot(3,1,2);
plot(EKF2.data(:,1)-t0,EKF2.data(:,4),'b');
grid on;
xlabel('time (sec)');ylabel('Y bias (cm/s^2)');
subplot(3,1,3);
plot(EKF2.data(:,1)-t0,EKF2.data(:,5),'b');
grid on;
xlabel('time (sec)');ylabel('Z bias (cm/s^2)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% Wind Velocity
fileName = 'WindVelEstimates';
figure;
subplot(2,1,1);
plot(EKF2.data(:,1)-t0,0.01*EKF2.data(:,6),'b');
grid on;
xlabel('time (sec)');ylabel('North Velocity (m/s)');
title('Wind Velocity Estimates');
subplot(2,1,2);
plot(EKF2.data(:,1)-t0,0.01*EKF2.data(:,7),'b');
grid on;
xlabel('time (sec)');ylabel('East Velocity (m/s)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% NED Magnetic Field
fileName = 'EarthMagEstimates';
figure;
subplot(3,1,1);
plot(EKF2.data(:,1)-t0,EKF2.data(:,8),'b');
grid on;
xlabel('time (sec)');ylabel('North Flux (mgauss)');
title('Earth Magnetic Field Estimates');
subplot(3,1,2);
plot(EKF2.data(:,1)-t0,EKF2.data(:,9),'b');
grid on;
xlabel('time (sec)');ylabel('East Flux (mgauss)');
subplot(3,1,3);
plot(EKF2.data(:,1)-t0,EKF2.data(:,10),'b');
grid on;
xlabel('time (sec)');ylabel('Down Flux (mgauss)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% XYZ Magnetic Field
fileName = 'BodyMagEstimates';
figure;
subplot(3,1,1);
plot(EKF2.data(:,1)-t0,EKF2.data(:,11),'b');
grid on;
xlabel('time (sec)');ylabel('X Flux (mgauss)');
title('Body Fixed Magnetic Field Bias Estimates');
subplot(3,1,2);
plot(EKF2.data(:,1)-t0,EKF2.data(:,12),'b');
grid on;
xlabel('time (sec)');ylabel('Y Flux (mgauss)');
subplot(3,1,3);
plot(EKF2.data(:,1)-t0,EKF2.data(:,13),'b');
grid on;
xlabel('time (sec)');ylabel('Z Flux (mgauss)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% Velocity Innovations
fileName = 'VelInnovations';
figure;
subplot(3,1,1);
plot(EKF3.data(:,1)-t0,0.01*[EKF3.data(:,3),-EKF4.data(:,3),EKF4.data(:,3)]);
grid on;
xlabel('time (sec)');ylabel('North (m/s)');
title('Velocity Measurement Innovations');
subplot(3,1,2);
plot(EKF3.data(:,1)-t0,0.01*[EKF3.data(:,4),-EKF4.data(:,4),EKF4.data(:,4)]);
grid on;
xlabel('time (sec)');ylabel('East (m/s)');
subplot(3,1,3);
plot(EKF3.data(:,1)-t0,0.01*[EKF3.data(:,5),-EKF4.data(:,5),EKF4.data(:,5)]);
grid on;
xlabel('time (sec)');ylabel('Down (m/s)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% Position Innovations
fileName = 'PosInnovations';
figure;
subplot(3,1,1);
plot(EKF3.data(:,1)-t0,0.01*[EKF3.data(:,6),-EKF4.data(:,6),EKF4.data(:,6)]);
grid on;
xlabel('time (sec)');ylabel('North (m)');
title('Position Measurement Innovations');
subplot(3,1,2);
plot(EKF3.data(:,1)-t0,0.01*[EKF3.data(:,7),-EKF4.data(:,7),EKF4.data(:,7)]);
grid on;
xlabel('time (sec)');ylabel('East (m)');
subplot(3,1,3);
plot(EKF3.data(:,1)-t0,0.01*[EKF3.data(:,8),-EKF4.data(:,8),EKF4.data(:,8)]);
grid on;
xlabel('time (sec)');ylabel('Down (m)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% Magnetometer Innovations
fileName = 'MagInnovations';
figure;
subplot(3,1,1);
plot(EKF3.data(:,1)-t0,[EKF3.data(:,9),-EKF4.data(:,9),EKF4.data(:,9)]);
grid on;
xlabel('time (sec)');ylabel('X Flux (mgauss)');
title('Magnetometer Measurement Innovations');
subplot(3,1,2);
plot(EKF3.data(:,1)-t0,[EKF3.data(:,10),-EKF4.data(:,10),EKF4.data(:,10)]);
grid on;
xlabel('time (sec)');ylabel('Y Flux (mgauss)');
subplot(3,1,3);
plot(EKF3.data(:,1)-t0,[EKF3.data(:,11),-EKF4.data(:,11),EKF4.data(:,11)]);
grid on;
xlabel('time (sec)');ylabel('Z Flux (mgauss)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% Airspeed Innovations
fileName = 'AirSpeedInnovations';
figure;
plot(EKF3.data(:,1)-t0,0.01*[EKF3.data(:,12),-EKF4.data(:,12),EKF4.data(:,12)]);
grid on;
xlabel('time (sec)');ylabel('airspeed (m/s)');
title('Airspeed Measurement Innovations');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');