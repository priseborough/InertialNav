%% plot gyro bias estimates
figure;
plot(statesLog(1,:),statesLog(8:10,:)/dt*180/pi);
grid on;
ylabel('Gyro Bias Estimate (deg/sec)');
xlabel('time (sec)');

%% plot velocity
figure;
plot(statesLog(1,:),statesLog(5:7,:));
grid on;
ylabel('Velocity (m/sec)');
xlabel('time (sec)');