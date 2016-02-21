function SaveScriptCodeTerrainEstimator3

%% Load Data
fileName = strcat('SymbolicOutput3.mat');
load(fileName);

%% Open output file
fileName = strcat('SymbolicOutput',int2str(nStates),'.txt');
fid = fopen(fileName,'wt');

%% Write equation for state transition matrix
if exist('F','var')
    
    fprintf(fid,'\n');
    fprintf(fid,'F = zeros(%d,%d);\n',nStates,nStates);
    for rowIndex = 1:nStates
        for colIndex = 1:nStates
            string = char(F(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'F(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end

%% Write equations for covariance prediction
% Only write out upper diagonal (matrix is symmetric)
if exist('SPP','var')
    
    fprintf(fid,'\n');
    fprintf(fid,'SPP = zeros(%d,1);\n',numel(SPP));
    for rowIndex = 1:numel(SPP)
        string = char(SPP(rowIndex,1));
        fprintf(fid,'SPP(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end

if exist('PP','var')
    
    fprintf(fid,'\n');
    fprintf(fid,'nextP = zeros(%d,%d);\n',nStates,nStates);
    for colIndex = 1:nStates
        for rowIndex = 1:colIndex
            string = char(PP(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'nextP(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end

%% Write equations for HAGL fusion
if exist('SK_HAGL','var')
    
    [nRow,nCol] = size(H_HAGL);
    fprintf(fid,'\n');
    fprintf(fid,'H_HAGL = zeros(1,%d);\n',nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(H_HAGL(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'H_HAGL(1,%d) = %s;\n',colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'\n');
    fprintf(fid,'SK_HAGL = zeros(%d,1);\n',numel(SK_HAGL));
    for rowIndex = 1:numel(SK_HAGL)
        string = char(SK_HAGL(rowIndex,1));
        fprintf(fid,'SK_HAGL(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
    [nRow,nCol] = size(K_HAGL);
    fprintf(fid,'\n');
    fprintf(fid,'Kfusion = zeros(%d,1);\n',nRow,nCol);
    for rowIndex = 1:nRow
        string = char(K_HAGL(rowIndex,1));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'Kfusion(%d) = %s;\n',rowIndex,string);
        end
    end
    fprintf(fid,'\n');
    
end
%% Write equations for vertical position fusion
if exist('SK_PD','var')
    
    [nRow,nCol] = size(H_HAGL);
    fprintf(fid,'\n');
    fprintf(fid,'H_PD = zeros(1,%d);\n',nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(H_PD(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'H_PD(1,%d) = %s;\n',colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'\n');
    fprintf(fid,'SK_PD = zeros(%d,1);\n',numel(SK_PD));
    for rowIndex = 1:numel(SK_PD)
        string = char(SK_PD(rowIndex,1));
        fprintf(fid,'SK_PD(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
    [nRow,nCol] = size(K_PD);
    fprintf(fid,'\n');
    fprintf(fid,'Kfusion = zeros(%d,1);\n',nRow,nCol);
    for rowIndex = 1:nRow
        string = char(K_PD(rowIndex,1));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'Kfusion(%d) = %s;\n',rowIndex,string);
        end
    end
    fprintf(fid,'\n');
    
end

%% Close output file
fclose(fid);

end