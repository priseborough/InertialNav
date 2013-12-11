maxTimeMag = IMUtime(find(IMUframe < max(MAGframe), 1, 'last' ));
maxTimeGps = IMUtime(find(IMUframe < max(GPSframe), 1, 'last' ));
maxTime = min(maxTimeMag,maxTimeGps);
dirName = '../plots/';

%% Euler Angles
euler_var = squeeze(euler_var);
fileName = 'EulerAngleEstimates';
figure;
subplot(3,1,1);
plot(time_est,[roll_est,euler_ref(:,1)]);
% hold on;
% plot(time_est,roll_est - sqrt(euler_var(:,1)),'r:');
% plot(time_est,roll_est + sqrt(euler_var(:,1)),'r:');
% hold off
grid on;
ylim([-200 200]);
xlabel('time (sec)');ylabel('roll (deg)');
title('Euler Angle Estimates');
subplot(3,1,2);
plot(time_est,[pitch_est,euler_ref(:,2)]);
% hold on;
% plot(time_est,pitch_est - sqrt(euler_var(:,2)),'r:');
% plot(time_est,pitch_est + sqrt(euler_var(:,2)),'r:');
% hold off
grid on;
ylim([-200 200]);
xlabel('time (sec)');ylabel('pitch (deg)');
subplot(3,1,3);
plot(time_est,[yaw_est,euler_ref(:,3)]);
% hold on;
% plot(time_est,yaw_est - sqrt(euler_var(:,3)),'r:');
% plot(time_est,yaw_est + sqrt(euler_var(:,3)),'r:');
% hold off
grid on;
ylim([-200 200]);
xlabel('time (sec)');ylabel('yaw (deg)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
saveas(gcf,strcat(dirName,fileName,'.pdf'));

%% NED velocity
fileName = 'VelocityEstimates';
figure;
subplot(3,1,1);
plot(time_est,[velN_est,velNED_ref(:,1)]);
hold on;
plot(time_est,velN_est+sqrt(variance_est(:,5)),'r:');
plot(time_est,velN_est-sqrt(variance_est(:,5)),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('North Velocity (m/s)');
title('NED Velocity Estimates');
subplot(3,1,2);
plot(time_est,[velE_est,velNED_ref(:,2)]);
hold on;
plot(time_est,velE_est+sqrt(variance_est(:,6)),'r:');
plot(time_est,velE_est-sqrt(variance_est(:,6)),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('East Velocity (m/s)');
subplot(3,1,3);
plot(time_est,[velD_est,velNED_ref(:,3)]);
hold on;
plot(time_est,velD_est+sqrt(variance_est(:,7)),'r:');
plot(time_est,velD_est-sqrt(variance_est(:,7)),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('Down Velocity (m/s)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
saveas(gcf,strcat(dirName,fileName,'.pdf'));

%% NE Position and Height
fileName = 'PositionEstimates';
figure;
subplot(3,1,1);
plot(time_est,[posN_est,posNED_ref(:,1)]);
hold on;
plot(time_est,posN_est+sqrt(variance_est(:,8)),'r:');
plot(time_est,posN_est-sqrt(variance_est(:,8)),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('North Position (m)');
title('NED Position Estimates');
subplot(3,1,2);
plot(time_est,[posE_est,posNED_ref(:,2)]);
hold on;
plot(time_est,posE_est+sqrt(variance_est(:,9)),'r:');
plot(time_est,posE_est-sqrt(variance_est(:,9)),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('East Position (m)');
subplot(3,1,3);
plot(time_est,-[posD_est,posNED_ref(:,3)]);
hold on;
plot(time_est,-posD_est+sqrt(variance_est(:,10)),'r:');
plot(time_est,-posD_est-sqrt(variance_est(:,10)),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('Height (m)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
saveas(gcf,strcat(dirName,fileName,'.pdf'));

%% Gyro Bias
fileName = 'GyroBiasEstimates';
figure;
subplot(3,1,1);
plot(time_est,dAngBiasX_est*50*60);
hold on;
plot(time_est,(dAngBiasX_est+sqrt(variance_est(:,11)))*50*60,'r:');
plot(time_est,(dAngBiasX_est-sqrt(variance_est(:,11)))*50*60,'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('X bias (deg/min)');
title('Gyro Bias Error Estimates');
subplot(3,1,2);
plot(time_est,dAngBiasY_est*50*60);
hold on;
plot(time_est,(dAngBiasY_est+sqrt(variance_est(:,12)))*50*60,'r:');
plot(time_est,(dAngBiasY_est-sqrt(variance_est(:,12)))*50*60,'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('Y bias (deg/min)');
subplot(3,1,3);
plot(time_est,dAngBiasZ_est*50*60);
hold on;
plot(time_est,(dAngBiasZ_est+sqrt(variance_est(:,13)))*50*60,'r:');
plot(time_est,(dAngBiasZ_est-sqrt(variance_est(:,13)))*50*60,'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('Z bias (deg/min)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
saveas(gcf,strcat(dirName,fileName,'.pdf'));

%% Wind Velocity
fileName = 'WindVelEstimates';
figure;
subplot(2,1,1);
plot(time_est,velWN_est);
hold on;
plot(time_est,(velWN_est+sqrt(variance_est(:,14))),'r:');
plot(time_est,(velWN_est-sqrt(variance_est(:,14))),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('North Velocity (m/s)');
title('Wind Velocity Estimates');
subplot(2,1,2);
plot(time_est,velWE_est);
hold on;
plot(time_est,(velWE_est+sqrt(variance_est(:,15))),'r:');
plot(time_est,(velWE_est-sqrt(variance_est(:,15))),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('East Velocity (m/s)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
saveas(gcf,strcat(dirName,fileName,'.pdf'));

%% NED Magnetic Field
fileName = 'EarthMagEstimates';
figure;
subplot(3,1,1);
plot(time_est,magN_est);
hold on;
plot(time_est,(magN_est+sqrt(variance_est(:,16))),'r:');
plot(time_est,(magN_est-sqrt(variance_est(:,16))),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('North Flux (mgauss)');
title('Earth Magnetic Field Estimates');
subplot(3,1,2);
plot(time_est,magE_est);
hold on;
plot(time_est,(magE_est+sqrt(variance_est(:,17))),'r:');
plot(time_est,(magE_est-sqrt(variance_est(:,17))),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('East Flux (mgauss)');
subplot(3,1,3);
plot(time_est,magD_est);
hold on;
plot(time_est,(magD_est+sqrt(variance_est(:,18))),'r:');
plot(time_est,(magD_est-sqrt(variance_est(:,18))),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('Down Flux (mgauss)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
saveas(gcf,strcat(dirName,fileName,'.pdf'));

%% XYZ Magnetic Field
fileName = 'BodyMagEstimates';
figure;
subplot(3,1,1);
plot(time_est,magX_est);
hold on;
plot(time_est,(magX_est+sqrt(variance_est(:,19))),'r:');
plot(time_est,(magX_est-sqrt(variance_est(:,19))),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('X Flux (mgauss)');
title('Body Fixed Magnetic Field Bias Estimates');
subplot(3,1,2);
plot(time_est,magY_est);
hold on;
plot(time_est,(magY_est+sqrt(variance_est(:,20))),'r:');
plot(time_est,(magY_est-sqrt(variance_est(:,20))),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('Y Flux (mgauss)');
subplot(3,1,3);
plot(time_est,magZ_est);
hold on;
plot(time_est,(magZ_est+sqrt(variance_est(:,21))),'r:');
plot(time_est,(magZ_est-sqrt(variance_est(:,21))),'r:');
hold off
grid on;
xlabel('time (sec)');ylabel('Z Flux (mgauss)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
saveas(gcf,strcat(dirName,fileName,'.pdf'));

%% Position & Height Innovations
fileName = 'PosInnovations';
figure;
subplot(3,1,1);
plot(time1,innovPosN);
grid on;
xlabel('time (sec)');ylabel('North (m)');
title('Position Measurement Innovations');
subplot(3,1,2);
plot(time1,innovPosE);
grid on;
xlabel('time (sec)');ylabel('East (m)');
subplot(3,1,3);
plot(time4,innovPosD);
grid on;
xlabel('time (sec)');ylabel('Down (m)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
saveas(gcf,strcat(dirName,fileName,'.pdf'));

%% Velocity Innovations
fileName = 'VelInovations';
figure;
subplot(3,1,1);
plot(time1,innovVelN);
grid on;
xlabel('time (sec)');ylabel('North (m/s)');
title('Velocity Measurement Innovations');
subplot(3,1,2);
plot(time1,innovVelE);
grid on;
xlabel('time (sec)');ylabel('East (m/s)');
subplot(3,1,3);
plot(time1,innovVelD);
grid on;
xlabel('time (sec)');ylabel('Down (m/s)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
saveas(gcf,strcat(dirName,fileName,'.pdf'));

%% Magnetometer Innovations
nSamples = min([length(innovMagX),length(innovMagY),length(innovMagZ)]);
fileName = 'MagInnovations';
figure;
subplot(3,1,1);
plot(time2(1:nSamples),innovMagX(1:nSamples,:));
grid on;
xlabel('time (sec)');ylabel('X Flux (mgauss)');
title('Magnetometer Measurement Innovations');
subplot(3,1,2);
plot(time2(1:nSamples),innovMagY(1:nSamples,:));
grid on;
xlabel('time (sec)');ylabel('Y Flux (mgauss)');
subplot(3,1,3);
plot(time2(1:nSamples),innovMagZ(1:nSamples,:));
grid on;
xlabel('time (sec)');ylabel('Z Flux (mgauss)');
saveas(gcf,strcat(dirName,fileName,'.fig'));
saveas(gcf,strcat(dirName,fileName,'.pdf'));

%% Airspeed Innovations
fileName = 'AirSpeedInnovations';
figure;
plot(time3,innovVtas);
grid on;
xlabel('time (sec)');ylabel('airspeed (m/s)');
title('Airspeed Measurement Innovations');
saveas(gcf,strcat(dirName,fileName,'.fig'));
saveas(gcf,strcat(dirName,fileName,'.pdf'));