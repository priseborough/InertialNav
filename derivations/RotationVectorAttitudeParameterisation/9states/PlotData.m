%% calculate and plot tilt correction magnitude
temp = sqrt(angErrLog(2,:).^2 + angErrLog(3,:).^2)*180/pi;
figure;
plot(angErrLog(1,:),temp);
grid on;
ylabel('Tilt correction magnitude (deg)');
xlabel('time (sec)');

%% plot gyro bias estimates
figure;
plot(statesLog(1,:),statesLog(8:10,:)/dt*180/pi);
grid on;
ylabel('Gyro Bias Estimate (deg/sec)');
xlabel('time (sec)');