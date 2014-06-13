clear all;
LoadUbuntuFlashData;
for i = 1:2
    
    % Find best vel delay
    maxIndex = 26;
    velDelay = zeros(1,maxIndex);
    velVar = velDelay;
    for index = 1:maxIndex
        msecVelDelay = (index-1)*20;
        save('timing.txt','alignTime','startTime','endTime','msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay','EAS2TAS','-ascii');
        dos('NavEKF');
        importfile('VelPosFuse.txt');
        velDelay(index) = msecVelDelay;
        velVar(index) = mean([var(VelPosFuse(:,2)),var(VelPosFuse(:,4)),var(VelPosFuse(:,6))]);
    end
    msecVelDelay = velDelay(velVar == min(velVar))
    
    % Find best pos delay
    maxIndex = 26;
    posDelay = zeros(1,maxIndex);
    posVar = posDelay;
    for index = 1:maxIndex
        msecPosDelay = (index-1)*20;
        save('timing.txt','alignTime','startTime','endTime','msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay','EAS2TAS','-ascii');
        dos('NavEKF');
        importfile('VelPosFuse.txt');
        posDelay(index) = msecPosDelay;
        posVar(index) = mean([var(VelPosFuse(:,8)),var(VelPosFuse(:,10))]);
    end
    msecPosDelay = posDelay(posVar == min(posVar))
    
    % Find best hgt delay
    maxIndex = 26;
    hgtDelay = zeros(1,maxIndex);
    hgtVar = hgtDelay;
    for index = 1:maxIndex
        msecHgtDelay = (index-1)*20;
        save('timing.txt','alignTime','startTime','endTime','msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay','EAS2TAS','-ascii');
        dos('NavEKF');
        importfile('VelPosFuse.txt');
        hgtDelay(index) = msecHgtDelay;
        hgtVar(index) = mean(var(VelPosFuse(:,12)));
    end
    msecHgtDelay = hgtDelay(hgtVar == min(hgtVar))
    
    % Find best magnetometer delay
    maxIndex = 11;
    magDelay = zeros(1,maxIndex);
    magVar = magDelay;
    for index = 1:maxIndex
        msecMagDelay = (index-1)*20;
        save('timing.txt','alignTime','startTime','endTime','msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay','EAS2TAS','-ascii');
        dos('NavEKF');
        importfile('MagFuse.txt');
        magDelay(index) = msecMagDelay;
        magVar(index) = mean([var(MagFuse(:,2)),var(MagFuse(:,4)),var(MagFuse(:,6))]);
    end
    msecMagDelay = magDelay(magVar == min(magVar))
    
    % Find best airspeed delay
    maxIndex = 21;
    tasDelay = zeros(1,maxIndex);
    tasVar = tasDelay;
    for index = 1:maxIndex
        msecTasDelay = (index-1)*20;
        save('timing.txt','alignTime','startTime','endTime','msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay','EAS2TAS','-ascii');
        dos('NavEKF');
        importfile('TasFuse.txt');
        tasDelay(index) = msecTasDelay;
        tasVar(index) = mean(var(TasFuse(:,2)));
    end
    msecTasDelay = tasDelay(tasVar == min(tasVar))
    
end
save('timing.txt','alignTime','startTime','endTime','msecVelDelay','msecPosDelay','msecHgtDelay','msecMagDelay','msecTasDelay','EAS2TAS','-ascii');
    figure;plot(velDelay,velVar);
    figure;plot(posDelay,posVar);
    figure;plot(hgtDelay,hgtVar);
    figure;plot(magDelay,magVar);
    figure;plot(tasDelay,tasVar);

%dos('NavFilterEKF')
%PlotCcodeOutput;