importfile('StateDataOut.txt');
importfile('EulDataOut.txt');
importfile('CovDataOut.txt');
importfile('RefVelPosDataOut.txt');
importfile('VelPosFuse.txt');
importfile('MagFuse.txt');
importfile('TasFuse.txt');

xmin = max([StateDataOut(1,1),EulDataOut(1,1),CovDataOut(1,1),RefVelPosDataOut(1,1),VelPosFuse(1,1),MagFuse(1,1),TasFuse(1,1)]);
xmax = min([max(StateDataOut(:,1)),max(EulDataOut(:,1)),max(CovDataOut(:,1)),max(RefVelPosDataOut(:,1)),max(VelPosFuse(:,1)),max(MagFuse(:,1)),max(TasFuse(:,1))]);
rad2deg = 180/pi;
close all;
%% Euler Angles
for i = 1:length(EulDataOut)
    if EulDataOut(i,7) > pi
        EulDataOut(i,7) = EulDataOut(i,7) - 2*pi;
    end
end
figure;
subplot(3,1,1);
plot(EulDataOut(:,1),EulDataOut(:,2:3)*rad2deg);
xlim([xmin,xmax]);
grid on;
ylim([-200 200]);
xlabel('time (sec)');ylabel('roll (deg)');
title('Euler Angle Estimates');
subplot(3,1,2);
plot(EulDataOut(:,1),EulDataOut(:,4:5)*rad2deg);
xlim([xmin,xmax]);
grid on;
ylim([-200 200]);
xlabel('time (sec)');ylabel('pitch (deg)');
subplot(3,1,3);
plot(EulDataOut(:,1),EulDataOut(:,6:7)*rad2deg);
xlim([xmin,xmax]);
grid on;
ylim([-200 200]);
xlabel('time (sec)');ylabel('yaw (deg)');

%% NED velocity
% remove repeated GPS data points
lastRow = RefVelPosDataOut(1,:);
for i = 2:length(RefVelPosDataOut)
    sameRow = ((sum(RefVelPosDataOut(i,:) == lastRow)) == 6);
    if sameRow
        RefVelPosDataOut(i,:) = NaN(1,7);
    else
        lastRow = RefVelPosDataOut(i,:);
    end
end
figure;
subplot(3,1,1);
plot(RefVelPosDataOut(:,1),RefVelPosDataOut(:,2),'.r');
xlim([xmin,xmax]);
hold on;
plot(StateDataOut(:,1),StateDataOut(:,6),'b');
hold off;
grid on;
xlabel('time (sec)');ylabel('North Velocity (m/s)');
title('NED Velocity Estimates');
subplot(3,1,2);
plot(RefVelPosDataOut(:,1),RefVelPosDataOut(:,3),'.r');
xlim([xmin,xmax]);
hold on;
plot(StateDataOut(:,1),StateDataOut(:,7),'b');
hold off;
grid on;
xlabel('time (sec)');ylabel('East Velocity (m/s)');
subplot(3,1,3);
plot(RefVelPosDataOut(:,1),RefVelPosDataOut(:,4),'.r');
xlim([xmin,xmax]);
hold on;
plot(StateDataOut(:,1),StateDataOut(:,8),'b');
hold off;
grid on;
xlabel('time (sec)');ylabel('Down Velocity (m/s)');

%% NE Position and Height
figure;
subplot(3,1,1);
plot(RefVelPosDataOut(:,1),RefVelPosDataOut(:,5),'.r');
xlim([xmin,xmax]);
hold on;
plot(StateDataOut(:,1),StateDataOut(:,9),'b');
hold off;
grid on;
xlabel('time (sec)');ylabel('North Position (m)');
title('NED Position Estimates');
subplot(3,1,2);
plot(RefVelPosDataOut(:,1),RefVelPosDataOut(:,6),'.r');
xlim([xmin,xmax]);
hold on;
plot(StateDataOut(:,1),StateDataOut(:,10),'b');
hold off;
grid on;
xlabel('time (sec)');ylabel('East Position (m)');
subplot(3,1,3);
plot(RefVelPosDataOut(:,1),RefVelPosDataOut(:,7),'.r');
xlim([xmin,xmax]);
hold on;
plot(StateDataOut(:,1),-StateDataOut(:,11),'b');
hold off;
grid on;
xlabel('time (sec)');ylabel('Height (m)');

%% Gyro Bias
figure;
dt = median(diff(StateDataOut(:,1)));
subplot(3,1,1);
plot(StateDataOut(:,1),StateDataOut(:,12)*rad2deg*60/dt);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('X bias (deg/min)');
title('Gyro Bias Error Estimates');
subplot(3,1,2);
plot(StateDataOut(:,1),StateDataOut(:,13)*rad2deg*60/dt);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('Y bias (deg/min)');
subplot(3,1,3);
plot(StateDataOut(:,1),StateDataOut(:,14)*rad2deg*60/dt);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('Z bias (deg/min)');

%% Accelerometer Bias
figure;
subplot(3,1,1);
plot(StateDataOut(:,1),StateDataOut(:,15)/dt);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('X bias (m/s^2)');
title('Accelerometer Bias Error Estimates');
subplot(3,1,2);
plot(StateDataOut(:,1),StateDataOut(:,16)/dt);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('Y bias (m/s^2)');
subplot(3,1,3);
plot(StateDataOut(:,1),StateDataOut(:,17)/dt);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('Z bias (m/s^2)');

%% Wind Velocity
figure;
subplot(2,1,1);
plot(StateDataOut(:,1),StateDataOut(:,18));
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('North Velocity (m/s)');
title('Wind Velocity Estimates');
subplot(2,1,2);
xlim([xmin,xmax]);
plot(StateDataOut(:,1),StateDataOut(:,19));
grid on;
xlabel('time (sec)');ylabel('East Velocity (m/s)');

%% NED Magnetic Field
figure;
subplot(3,1,1);
plot(StateDataOut(:,1),StateDataOut(:,20)*1024);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('North Flux (mgauss)');
title('Earth Magnetic Field Estimates');
subplot(3,1,2);
plot(StateDataOut(:,1),StateDataOut(:,21)*1024);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('East Flux (mgauss)');
subplot(3,1,3);
plot(StateDataOut(:,1),StateDataOut(:,22)*1024);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('Down Flux (mgauss)');

%% XYZ Magnetic Field
figure;
subplot(3,1,1);
plot(StateDataOut(:,1),StateDataOut(:,23)*1024);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('X Flux (mgauss)');
title('Body Fixed Magnetic Field Bias Estimates');
subplot(3,1,2);
plot(StateDataOut(:,1),StateDataOut(:,24)*1024);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('Y Flux (mgauss)');
subplot(3,1,3);
plot(StateDataOut(:,1),StateDataOut(:,25)*1024);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('Z Flux (mgauss)');

%% Velocity Innovations
figure;
subplot(3,1,1);
plot(VelPosFuse(:,1),[VelPosFuse(:,2),-sqrt(VelPosFuse(:,3)),sqrt(VelPosFuse(:,3))]);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('North (m/s)');
title('Velocity Measurement Innovations');
subplot(3,1,2);
plot(VelPosFuse(:,1),[VelPosFuse(:,4),-sqrt(VelPosFuse(:,5)),sqrt(VelPosFuse(:,5))]);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('East (m/s)');
subplot(3,1,3);
plot(VelPosFuse(:,1),[VelPosFuse(:,6),-sqrt(VelPosFuse(:,7)),sqrt(VelPosFuse(:,7))]);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('Down (m/s)');

%% Position Innovations
figure;
subplot(3,1,1);
plot(VelPosFuse(:,1),[VelPosFuse(:,8),-sqrt(VelPosFuse(:,9)),sqrt(VelPosFuse(:,9))]);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('North (m)');
title('Position Measurement Innovations');
subplot(3,1,2);
plot(VelPosFuse(:,1),[VelPosFuse(:,10),-sqrt(VelPosFuse(:,11)),sqrt(VelPosFuse(:,11))]);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('East (m)');
subplot(3,1,3);
plot(VelPosFuse(:,1),[VelPosFuse(:,12),-sqrt(VelPosFuse(:,13)),sqrt(VelPosFuse(:,13))]);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('Down (m)');

%% Magnetometer Innovations
figure;
subplot(3,1,1);
plot(MagFuse(:,1),[MagFuse(:,2),-sqrt(MagFuse(:,3)),sqrt(MagFuse(:,3))]*1024);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('X Flux (mgauss)');
title('Magnetometer Measurement Innovations');
subplot(3,1,2);
plot(MagFuse(:,1),[MagFuse(:,4),-sqrt(MagFuse(:,5)),sqrt(MagFuse(:,5))]*1024);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('Y Flux (mgauss)');
subplot(3,1,3);
plot(MagFuse(:,1),[MagFuse(:,6),-sqrt(MagFuse(:,7)),sqrt(MagFuse(:,7))]*1024);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('Z Flux (mgauss)');

%% AHRS Diagnostic
xmin = 310;
xmax = 330;
figure;
subplot(4,1,1);
plot(EulDataOut(:,1),EulDataOut(:,4:5)*rad2deg);
xlim([xmin,xmax]);
grid on;
xlabel('time (sec)');ylabel('pitch (deg)');
title('AHRS Diagnostic');
subplot(4,1,2);
plot(RefVelPosDataOut(:,1),RefVelPosDataOut(:,7),'r');
xlim([xmin,xmax]);
hold on;
plot(StateDataOut(:,1),-StateDataOut(:,11),'b');
hold off;
grid on;
xlabel('time (sec)');ylabel('Height (m)');
subplot(4,1,3);
plot(EulDataOut(:,1),EulDataOut(:,8));
xlim([xmin,xmax]);
grid on;
ylim([0 1.0]);
xlabel('time (sec)');ylabel('AHRS Error RP');
subplot(4,1,4);
plot(EulDataOut(:,1),EulDataOut(:,9));
xlim([xmin,xmax]);
grid on;
ylim([0 1.0]);
xlabel('time (sec)');ylabel('AHRS Error Yaw');
