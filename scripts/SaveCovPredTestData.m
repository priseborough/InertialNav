Pin=double(Pin);
save('../code/Pin.txt','Pin','-ascii','-double');
Pout=double(Pout);
save('../code/Pout.txt','Pout','-ascii','-double');
correctedDelAng=double(correctedDelAng);
save('../code/correctedDelAng.txt','correctedDelAng','-ascii','-double');
correctedDelVel=double(correctedDelVel);
save('../code/correctedDelVel.txt','correctedDelVel','-ascii','-double');
states=double(states);
save('../code/states.txt','states','-ascii','-double');
dt=double(dt);
save('../code/dt.txt','dt','-ascii','-double');
onGround=double(onGround);
save('../code/onGround.txt','onGround','-ascii','-double');
useAirspeed=double(useAirspeed);
save('../code/useAirspeed.txt','useAirspeed','-ascii','-double');
useCompass=double(useCompass);
save('../code/useCompass.txt','useCompass','-ascii','-double');

save('CovPredTestVectors.mat','correctedDelAng','correctedDelVel','dt','Pin','Pout','states','onGround','useAirspeed','useCompass');

