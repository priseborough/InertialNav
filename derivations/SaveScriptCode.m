function SaveScriptCode(nStates)
%% Load Data
fileName = strcat('SymbolicOutput',int2str(nStates),'.mat');
load(fileName);

%% Open output file
fileName = strcat('SymbolicOutput',int2str(nStates),'.txt');
fid = fopen(fileName,'wt');

%% Write equation for state transition matrix
fprintf(fid,'SF = zeros(%d,1);\n',numel(SF));
for rowIndex = 1:numel(SF)
    string = char(SF(rowIndex,1));
    fprintf(fid,'SF(%d) = %s;\n',rowIndex,string);
end

% fprintf(fid,'\n');
% fprintf(fid,'F = zeros(%d,%d);\n',nStates,nStates);
% for rowIndex = 1:nStates
%     for colIndex = 1:nStates
%         string = char(F(rowIndex,colIndex));
%         % don't write out a zero-assignment
%         if ~strcmpi(string,'0')
%             fprintf(fid,'F(%d,%d) = %s;\n',rowIndex,colIndex,string);
%         end
%     end
% end
% fprintf(fid,'\n');

%% Write equations for control influence (disturbance) matrix
fprintf(fid,'\n');
fprintf(fid,'SG = zeros(%d,1);\n',numel(SG));
for rowIndex = 1:numel(SG)
    string = char(SG(rowIndex,1));
    fprintf(fid,'SG(%d) = %s;\n',rowIndex,string);
end
fprintf(fid,'\n');

% fprintf(fid,'\n');
% fprintf(fid,'G = zeros(%d,%d);\n',nStates,numel([da;dv]));
% for rowIndex = 1:nStates
%     for colIndex = 1:numel([da;dv])
%         string = char(G(rowIndex,colIndex));
%         % don't write out a zero-assignment
%         if ~strcmpi(string,'0')
%             fprintf(fid,'G(%d,%d) = %s;\n',rowIndex,colIndex,string);
%         end
%     end
% end
% fprintf(fid,'\n');

%% Write equations for state error matrix
fprintf(fid,'\n');
fprintf(fid,'SQ = zeros(%d,1);\n',numel(SQ));
for rowIndex = 1:numel(SQ)
    string = char(SQ(rowIndex,1));
    fprintf(fid,'SQ(%d) = %s;\n',rowIndex,string);
end
fprintf(fid,'\n');

% fprintf(fid,'\n');
% fprintf(fid,'Q = zeros(%d,%d);\n',nStates,nStates);
% for rowIndex = 1:nStates
%     for colIndex = 1:nStates
%         string = char(Q(rowIndex,colIndex));
%         % don't write out a zero-assignment
%         if ~strcmpi(string,'0')
%             fprintf(fid,'Q(%d,%d) = %s;\n',rowIndex,colIndex,string);
%         end
%     end
% end
% fprintf(fid,'\n');

%% Write equations for covariance prediction
fprintf(fid,'\n');
fprintf(fid,'SPP = zeros(%d,1);\n',numel(SPP));
for rowIndex = 1:numel(SPP)
    string = char(SPP(rowIndex,1));
    fprintf(fid,'SPP(%d) = %s;\n',rowIndex,string);
end
fprintf(fid,'\n');

fprintf(fid,'\n');
fprintf(fid,'nextP = zeros(%d,%d);\n',nStates,nStates);
for rowIndex = 1:nStates
    for colIndex = 1:nStates
        string = char(PP(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'nextP(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end
fprintf(fid,'\n');

%% Write equations for velocity and position data fusion
[nRow,nCol] = size(H_VP);
fprintf(fid,'\n');
fprintf(fid,'H_VP = zeros(%d,%d);\n',nRow,nCol);
for rowIndex = 1:nRow
    for colIndex = 1:nCol
        string = char(H_VP(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'H_VP(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end
fprintf(fid,'\n');

[nRow,nCol] = size(SK_VP);
fprintf(fid,'\n');
fprintf(fid,'SK_VP = zeros(%d,%d);\n',nRow,nCol);
for rowIndex = 1:nRow
    for colIndex = 1:nCol
        string = char(SK_VP(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'SK_VP(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end
fprintf(fid,'\n');

[nRow,nCol] = size(K_VP);
fprintf(fid,'\n');
fprintf(fid,'Kfusion = zeros(%d,%d);\n',nRow,nCol);
for rowIndex = 1:nRow
    for colIndex = 1:nCol
        string = char(K_VP(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'Kfusion(%d,%d) = %s;\n',rowIndex,colIndex,string);
        end
    end
end
fprintf(fid,'\n');

%% Write equations for true airspeed data fusion
fprintf(fid,'\n');
fprintf(fid,'SH_TAS = zeros(%d,1);\n',numel(SH_TAS));
for rowIndex = 1:numel(SH_TAS)
    string = char(SH_TAS(rowIndex,1));
    fprintf(fid,'SH_TAS(%d) = %s;\n',rowIndex,string);
end

[nRow,nCol] = size(H_TAS);
fprintf(fid,'\n');
fprintf(fid,'H_TAS = zeros(1,%d);\n',nCol);
for rowIndex = 1:nRow
    for colIndex = 1:nCol
        string = char(H_TAS(rowIndex,colIndex));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'H_TAS(1,%d) = %s;\n',colIndex,string);
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'\n');
fprintf(fid,'SK_TAS = zeros(%d,1);\n',numel(SK_TAS));
for rowIndex = 1:numel(SK_TAS)
    string = char(SK_TAS(rowIndex,1));
    fprintf(fid,'SK_TAS(%d) = %s;\n',rowIndex,string);
end
fprintf(fid,'\n');

[nRow,nCol] = size(K_TAS);
fprintf(fid,'\n');
fprintf(fid,'Kfusion = zeros(%d,1);\n',nRow,nCol);
for rowIndex = 1:nRow
    string = char(K_TAS(rowIndex,1));
    % don't write out a zero-assignment
    if ~strcmpi(string,'0')
        fprintf(fid,'Kfusion(%d) = %s;\n',rowIndex,string);
    end
end
fprintf(fid,'\n');

%% Write equations for magnetometer data fusion
fprintf(fid,'\n');
fprintf(fid,'SH_MAG = zeros(%d,1);\n',numel(SH_MAG));
for rowIndex = 1:numel(SH_MAG)
    string = char(SH_MAG(rowIndex,1));
    fprintf(fid,'SH_MAG(%d) = %s;\n',rowIndex,string);
end
fprintf(fid,'\n');

[nRow,nCol] = size(H_MAG);
fprintf(fid,'\n');
fprintf(fid,'H_MAG = zeros(1,%d);\n',nCol);
for colIndex = 1:nCol
    string = char(H_MAG(1,colIndex));
    % don't write out a zero-assignment
    if ~strcmpi(string,'0')
        fprintf(fid,'H_MAG(%d) = %s;\n',colIndex,string);
    end
end
fprintf(fid,'\n');

fprintf(fid,'\n');
fprintf(fid,'SK_MX = zeros(%d,1);\n',numel(SK_MX));
for rowIndex = 1:numel(SK_MX)
    string = char(SK_MX(rowIndex,1));
    fprintf(fid,'SK_MX(%d) = %s;\n',rowIndex,string);
end
fprintf(fid,'\n');

[nRow,nCol] = size(K_MX);
fprintf(fid,'\n');
fprintf(fid,'Kfusion = zeros(%d,1);\n',nRow,nCol);
for rowIndex = 1:nRow
    string = char(K_MX(rowIndex,1));
    % don't write out a zero-assignment
    if ~strcmpi(string,'0')
        fprintf(fid,'Kfusion(%d) = %s;\n',rowIndex,string);
    end
end
fprintf(fid,'\n');

[nRow,nCol] = size(H_MAG);
fprintf(fid,'\n');
fprintf(fid,'H_MAG = zeros(1,%d);\n',nCol);
for colIndex = 1:nCol
    string = char(H_MAG(2,colIndex));
    % don't write out a zero-assignment
    if ~strcmpi(string,'0')
        fprintf(fid,'H_MAG(%d) = %s;\n',colIndex,string);
    end
end
fprintf(fid,'\n');

fprintf(fid,'\n');
fprintf(fid,'SK_MY = zeros(%d,1);\n',numel(SK_MY));
for rowIndex = 1:numel(SK_MY)
    string = char(SK_MY(rowIndex,1));
    fprintf(fid,'SK_MY(%d) = %s;\n',rowIndex,string);
end
fprintf(fid,'\n');

[nRow,nCol] = size(K_MY);
fprintf(fid,'\n');
fprintf(fid,'Kfusion = zeros(%d,1);\n',nRow,nCol);
for rowIndex = 1:nRow
    string = char(K_MY(rowIndex,1));
    % don't write out a zero-assignment
    if ~strcmpi(string,'0')
        fprintf(fid,'Kfusion(%d) = %s;\n',rowIndex,string);
    end
end
fprintf(fid,'\n');

[nRow,nCol] = size(H_MAG);
fprintf(fid,'\n');
fprintf(fid,'H_MAG = zeros(1,%d);\n',nCol);
for colIndex = 1:nCol
    string = char(H_MAG(3,colIndex));
    % don't write out a zero-assignment
    if ~strcmpi(string,'0')
        fprintf(fid,'H_MAG(%d) = %s;\n',colIndex,string);
    end
end
fprintf(fid,'\n');

fprintf(fid,'\n');
fprintf(fid,'SK_MZ = zeros(%d,1);\n',numel(SK_MZ));
for rowIndex = 1:numel(SK_MZ)
    string = char(SK_MZ(rowIndex,1));
    fprintf(fid,'SK_MZ(%d) = %s;\n',rowIndex,string);
end
fprintf(fid,'\n');

[nRow,nCol] = size(K_MZ);
fprintf(fid,'\n');
fprintf(fid,'Kfusion = zeros(%d,1);\n',nRow,nCol);
for rowIndex = 1:nRow
    string = char(K_MZ(rowIndex,1));
    % don't write out a zero-assignment
    if ~strcmpi(string,'0')
        fprintf(fid,'Kfusion(%d) = %s;\n',rowIndex,string);
    end
end
fprintf(fid,'\n');

%% Write equations for optical flow sensor angular LOS data fusion
fprintf(fid,'\n');
fprintf(fid,'SH_LOS = zeros(%d,1);\n',numel(SH_LOS));
for rowIndex = 1:numel(SH_LOS)
    string = char(SH_LOS(rowIndex,1));
    fprintf(fid,'SH_LOS(%d) = %s;\n',rowIndex,string);
end
fprintf(fid,'\n');

fprintf(fid,'\n');
[nRow,nCol] = size(H_LOS);
fprintf(fid,'\n');
fprintf(fid,'H_LOSX = zeros(1,%d);\n',nCol);
for colIndex = 1:nCol
    string = char(H_LOS(1,colIndex));
    % don't write out a zero-assignment
    if ~strcmpi(string,'0')
        fprintf(fid,'H_LOSX(%d) = %s;\n',colIndex,string);
    end
end
fprintf(fid,'\n');

fprintf(fid,'\n');
[nRow,nCol] = size(H_LOS);
fprintf(fid,'\n');
fprintf(fid,'H_LOSY = zeros(1,%d);\n',nCol);
for colIndex = 1:nCol
    string = char(H_LOS(2,colIndex));
    % don't write out a zero-assignment
    if ~strcmpi(string,'0')
        fprintf(fid,'H_LOSY(%d) = %s;\n',colIndex,string);
    end
end

fprintf(fid,'\n');
fprintf(fid,'SK_LOS = zeros(%d,1);\n',numel(SK_LOS));
for rowIndex = 1:numel(SK_LOS)
    string = char(SK_LOS(rowIndex,1));
    fprintf(fid,'SK_LOS(%d) = %s;\n',rowIndex,string);
end

[nRow,nCol] = size(K_LOSX);
fprintf(fid,'\n');
fprintf(fid,'Kfusion = zeros(%d,1);\n',nRow,nCol);
for rowIndex = 1:nRow
    string = char(K_LOS(rowIndex,1));
    % don't write out a zero-assignment
    if ~strcmpi(string,'0')
        fprintf(fid,'Kfusion(%d) = %s;\n',rowIndex,string);
    end
end
fprintf(fid,'\n');

[nRow,nCol] = size(K_LOSY);
fprintf(fid,'\n');
fprintf(fid,'Kfusion = zeros(%d,1);\n',nRow,nCol);
for rowIndex = 1:nRow
    string = char(K_LOS(rowIndex,2));
    % don't write out a zero-assignment
    if ~strcmpi(string,'0')
        fprintf(fid,'Kfusion(%d) = %s;\n',rowIndex,string);
    end
end

%% Close output file
fclose(fid);

end