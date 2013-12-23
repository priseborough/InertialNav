fid = fopen('Pin.txt', 'r');
dataIn = fscanf (fid, '%lf');
fclose(fid);
Pin = zeros(24,24);
index = 1;
for i=1:24
    for j=1:24
        Pin(i,j) = dataIn(index);
        index = index + 1;
    end
end

fid = fopen('Pout.txt', 'r');
dataIn = fscanf (fid, '%lf');
fclose(fid);
Pref = zeros(24,24);
index = 1;
for i=1:24
    for j=1:24
        Pref(i,j) = dataIn(index);
        index = index + 1;
    end
end

fid = fopen('states.txt', 'r');
dataIn = fscanf (fid, '%lf');
fclose(fid);
for i=1:24
    states(i) = dataIn(i);
end

fid = fopen('correctedDelAng.txt', 'r');
dataIn = fscanf (fid, '%lf');
fclose(fid);
for i=1:3
    correctedDelAng(i) = dataIn(i);
end

fid = fopen('dt.txt', 'r');
dt = fscanf (fid, '%lf');
fclose(fid);

fid = fopen('onGround.txt', 'r');
onGround = boolean(fscanf (fid, '%lf'));
fclose(fid);

fid = fopen('useAirspeed.txt', 'r');
useAirspeed = boolean(fscanf (fid, '%lf'));
fclose(fid);

fid = fopen('useCompass.txt', 'r');
useCompass = boolean(fscanf (fid, '%lf'));
fclose(fid);