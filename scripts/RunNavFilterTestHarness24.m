clear all;
close all;
LoadNavFilterTestDataStruct
alignTime = min(IMUtime(IMUframe>GPSframe(find(GndSpd  >8, 1 )))) - 10;
startTime = alignTime - 30;
endTime = max(IMUtime)-1;
msecVelDelay = 200;
msecPosDelay = 160;
msecHgtDelay = 340;
msecMagDelay = 20;
msecTasDelay = 220;
EAS2TAS = 1.08;
sim('NavFilterTestHarness24')
PlotNavFilterData24