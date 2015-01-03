function SaveScriptCode
%% Load Data
nStates = 1;
fileName = strcat('SymbolicOutput.mat');
load(fileName);

%% Open output file
fileName = strcat('SymbolicOutput.txt');
fid = fopen(fileName,'wt');

%% Write equation for state transition matrix
if exist('SF','var')
    
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
    
end
%% Write equations for control influence (disturbance) matrix
if exist('SG','var')
    
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
    
end
%% Write equations for state error matrix
if exist('SQ','var')
    
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
    
end
%% Write equations for covariance prediction
if exist('nextPopt','var')
    
    fprintf(fid,'\n');
    fprintf(fid,'nextPopt = zeros(%d,%d);\n',nStates,nStates);
    for rowIndex = 1:nStates
        for colIndex = 1:nStates
            string = char(nextPopt(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'nextPopt(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end

%% Write equations for optical flow rate fusion

% Kalman Gain
if exist('SH_OPT','var')
    
    fprintf(fid,'\n');
    fprintf(fid,'SH_OPT = zeros(%d,1);\n',numel(SH_OPT));
    for rowIndex = 1:numel(SH_OPT)
        string = char(SH_OPT(rowIndex,1));
        fprintf(fid,'SH_OPT(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'\n');
    [nRow,nCol] = size(H_OPT);
    fprintf(fid,'H_OPT = zeros(%d,%d);\n',nRow,nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(H_OPT(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'H_OPT(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');

    fprintf(fid,'\n');
    fprintf(fid,'SK_OPT = zeros(%d,1);\n',numel(SK_OPT));
    for rowIndex = 1:numel(SK_OPT)
        string = char(SK_OPT(rowIndex,1));
        fprintf(fid,'SK_OPT(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'\n');
    [nRow,nCol] = size(K_OPT);
    fprintf(fid,'K_OPT = zeros(%d,%d);\n',nRow,nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(K_OPT(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'K_OPT(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end

% X rate covariance update
if exist('SPX','var')
    
    fprintf(fid,'\n');
    fprintf(fid,'SPX = zeros(%d,1);\n',numel(SPX));
    for rowIndex = 1:numel(SPX)
        string = char(SPX(rowIndex,1));
        fprintf(fid,'SPX(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end

if exist('nextPX','var')
    
    fprintf(fid,'\n');
    fprintf(fid,'nextPopt = zeros(%d,%d);\n',nStates,nStates);
    for rowIndex = 1:nStates
        for colIndex = 1:nStates
            string = char(nextPX(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'nextPopt(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end

% Y rate covariance update
if exist('SPY','var')
    
    fprintf(fid,'\n');
    fprintf(fid,'SPY = zeros(%d,1);\n',numel(SPY));
    for rowIndex = 1:numel(SPY)
        string = char(SPY(rowIndex,1));
        fprintf(fid,'SPY(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end

if exist('nextPY','var')
    
    fprintf(fid,'\n');
    fprintf(fid,'nextPopt = zeros(%d,%d);\n',nStates,nStates);
    for rowIndex = 1:nStates
        for colIndex = 1:nStates
            string = char(nextPY(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'nextPopt(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end

%% Write equations for range fusion
if exist('H_RNG','var')
    
    [nRow,nCol] = size(H_RNG);
    fprintf(fid,'\n');
    fprintf(fid,'H_RNG = zeros(1,%d);\n',nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(H_RNG(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'H_RNG(1,%d) = %s;\n',colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'\n');
    fprintf(fid,'SK_RNG = zeros(%d,1);\n',numel(SK_RNG));
    for rowIndex = 1:numel(SK_RNG)
        string = char(SK_RNG(rowIndex,1));
        fprintf(fid,'SK_RNG(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
    [nRow,nCol] = size(K_RNG);
    fprintf(fid,'\n');
    fprintf(fid,'K_RNG = zeros(%d,1);\n',nRow,nCol);
    for rowIndex = 1:nRow
        string = char(K_RNG(rowIndex,1));
        % don't write out a zero-assignment
        if ~strcmpi(string,'0')
            fprintf(fid,'K_RNG(%d) = %s;\n',rowIndex,string);
        end
    end
    fprintf(fid,'\n');
    
end

if exist('SPR','var')
    
    fprintf(fid,'\n');
    fprintf(fid,'SPR = zeros(%d,1);\n',numel(SPR));
    for rowIndex = 1:numel(SPR)
        string = char(SPR(rowIndex,1));
        fprintf(fid,'SPR(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end

if exist('nextPR','var')
    
    fprintf(fid,'\n');
    fprintf(fid,'nextPopt = zeros(%d,%d);\n',nStates,nStates);
    for rowIndex = 1:nStates
        for colIndex = 1:nStates
            string = char(nextPR(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'nextPopt(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end

%% Close output file
fclose(fid);

end
