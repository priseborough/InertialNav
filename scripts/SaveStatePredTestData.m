statesInINS=double(statesInINS);
save('../code/statesInINS.txt','statesInINS','-ascii','-double');
earthRateNED=double(earthRateNED);
save('../code/earthRateNED.txt','earthRateNED','-ascii','-double');
dAngIMU=double(dAngIMU);
save('../code/dAngIMU.txt','dAngIMU','-ascii','-double');
dVelIMU=double(dVelIMU);
save('../code/dVelIMU.txt','dVelIMU','-ascii','-double');
dtIMU=double(dtIMU);
save('../code/dtIMU.txt','dtIMU','-ascii','-double');
statesOutINS=double(statesOutINS);
save('../code/statesOutINS.txt','statesOutINS','-ascii','-double');

save('StatePredTestVectors.mat','statesInINS','earthRateNED','dAngIMU','dVelIMU','dtIMU','statesOutINS');

