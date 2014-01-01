%% Load data
load('simtest9.mat');
dirName = '../plots/';
t0 = EKF1.data(1,1);
%% Euler Angles
fileName = 'EulerAngleEstimates';
figure;
subplot(3,1,1);
plot(ATT.data(:,1)-t0,ATT.data(:,3),'g');
hold on;
plot(EKF1.data(:,1)-t0,EKF1.data(:,3),'b');
hold off;
grid on;
ylim([-180 180]);
xlabel('time (sec)');ylabel('roll (deg)');
title('Euler Angle Estimates');
subplot(3,1,2);
plot(ATT.data(:,1)-t0,ATT.data(:,4),'g');
hold on;
plot(EKF1.data(:,1)-t0,EKF1.data(:,4),'b');
hold off;
grid on;
ylim([-180 180]);
xlabel('time (sec)');ylabel('pitch (deg)');
subplot(3,1,3);
plot(ATT.data(:,1)-t0,ATT.data(:,5),'g');
hold on;
plot(EKF1.data(:,1)-t0,EKF1.data(:,5),'b');
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
plot(GPS.data(:,1)-t0,GPS.data(:,11).*cos(GPS.data(:,12)*pi/180),'g');
hold on;
plot(EKF1.data(:,1)-t0,EKF1.data(:,6),'b');
hold off;
grid on;
xlabel('time (sec)');ylabel('North Velocity (m/s)');
title('NED Velocity Estimates');
subplot(3,1,2);
plot(GPS.data(:,1)-t0,GPS.data(:,11).*sin(GPS.data(:,12)*pi/180),'g');
hold on;
plot(EKF1.data(:,1)-t0,EKF1.data(:,7),'b');
hold off;
hold off;
grid on;
xlabel('time (sec)');ylabel('East Velocity (m/s)');
subplot(3,1,3);
plot(GPS.data(:,1)-t0,GPS.data(:,13),'g');
hold on;
plot(EKF1.data(:,1)-t0,EKF1.data(:,8),'b');
hold off;
grid on;
xlabel('time (sec)');ylabel('Down Velocity (m/s)');
saveas(gcf,strcat(dirName,fileName,'.fig')); 
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% NE Position and Height
fileName = 'PositionEstimates';
figure;
subplot(3,1,1);
plot(EKF1.data(:,1)-t0,EKF1.data(:,9),'b');
grid on;
xlabel('time (sec)');ylabel('North Position (m)');
title('NED Position Estimates');
subplot(3,1,2);
plot(EKF1.data(:,1)-t0,EKF1.data(:,10),'b');
grid on;
xlabel('time (sec)');ylabel('East Position (m)');
subplot(3,1,3);
plot(GPS.data(:,1)-t0,GPS.data(:,9)-GPS.data(find(GPS.data(:,9)~=0, 1 ),9),'g');
hold on;
plot(EKF1.data(:,1)-t0,-EKF1.data(:,11),'b');
hold off
grid on;
ylim([-10,max(-EKF1.data(:,11))+10]);
xlabel('time (sec)');ylabel('Height (m)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% Gyro Bias
fileName = 'GyroBiasEstimates';
figure;
subplot(3,1,1);
plot(EKF1.data(:,1)-t0,EKF1.data(:,12),'b');
grid on;
xlabel('time (sec)');ylabel('X bias (deg/min)');
title('Gyro Bias Error Estimates');
subplot(3,1,2);
plot(EKF1.data(:,1)-t0,EKF1.data(:,13),'b');
grid on;
xlabel('time (sec)');ylabel('Y bias (deg/min)');
subplot(3,1,3);
plot(EKF1.data(:,1)-t0,EKF1.data(:,14),'b');
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
xlabel('time (sec)');ylabel('X bias (m/s^2)');
title('Accelerometer Bias Error Estimates');
subplot(3,1,2);
plot(EKF2.data(:,1)-t0,EKF2.data(:,4),'b');
grid on;
xlabel('time (sec)');ylabel('Y bias (m/s^2)');
subplot(3,1,3);
plot(EKF2.data(:,1)-t0,EKF2.data(:,5),'b');
grid on;
xlabel('time (sec)');ylabel('Z bias (m/s^2)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
print(gcf, '-djpeg', strcat(dirName,fileName,'.jpg'), '-r200');

%% Wind Velocity
fileName = 'WindVelEstimates';
figure;
subplot(2,1,1);
plot(EKF2.data(:,1)-t0,EKF2.data(:,6),'b');
grid on;
xlabel('time (sec)');ylabel('North Velocity (m/s)');
title('Wind Velocity Estimates');
subplot(2,1,2);
plot(EKF2.data(:,1)-t0,EKF2.data(:,7),'b');
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