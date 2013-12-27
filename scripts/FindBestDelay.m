clear all;
LoadNavFilterTestDataStruct;
alignTime = min(IMUtime(IMUframe>GPSframe(find(GndSpd  >8, 1 )))) - 10;
startTime = alignTime - 30;
endTime = max(IMUtime)-1;
msecVelDelay = 180;
msecPosDelay = 160;
msecHgtDelay = 360;
msecMagDelay = 20;
msecTasDelay = 280;
EAS2TAS = 1.08;
for i = 1:1
    
    % Find best vel delay
    maxIndex = 21;
    delay = zeros(1,maxIndex);
    variance = delay;
    for index = 1:maxIndex
        msecVelDelay = (index-1)*20;
        sim('NavFilterTestHarness24')
        delay(index) = msecVelDelay;
        variance(index) = mean([var(innovVelN),var(innovVelE),var(innovVelD)]);
    end
    figure;plot(delay,variance);
    msecVelDelay = delay(variance == min(variance))
    
    % Find best pos delay
    maxIndex = 21;
    delay = zeros(1,maxIndex);
    variance = delay;
    for index = 1:maxIndex
        msecPosDelay = (index-1)*20;
        sim('NavFilterTestHarness24')
        delay(index) = msecPosDelay;
        variance(index) = mean([var(innovPosN),var(innovPosE)]);
    end
    figure;plot(delay,variance);
    msecPosDelay = delay(variance == min(variance))
    % Find best hgt delay
    maxIndex = 26;
    delay = zeros(1,maxIndex);
    variance = delay;
    for index = 1:maxIndex
        msecHgtDelay = (index-1)*20;
        sim('NavFilterTestHarness24')
        delay(index) = msecHgtDelay;
        variance(index) = mean(var(innovPosD));
    end
    figure;plot(delay,variance);
    msecHgtDelay = delay(variance == min(variance))
    
    % Find best magnetometer delay
    maxIndex = 11;
    delay = zeros(1,maxIndex);
    variance = delay;
    for index = 1:maxIndex
        msecMagDelay = (index-1)*20;
        sim('NavFilterTestHarness24')
        delay(index) = msecMagDelay;
        variance(index) = mean([var(innovMagX),var(innovMagY),var(innovMagZ)]);
    end
    figure;plot(delay,variance);
    msecMagDelay = delay(variance == min(variance))
    
    % Find best airspeed delay
    maxIndex = 20;
    delay = zeros(1,maxIndex);
    variance = delay;
    for index = 1:maxIndex
        msecTasDelay = (index-1)*20;
        sim('NavFilterTestHarness24')
        delay(index) = msecTasDelay;
        variance(index) = mean(var(innovVtas));
    end
    figure;plot(delay,variance);
    msecTasDelay = delay(variance == min(variance))
    
end
PlotNavFilter24