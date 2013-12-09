%% Initialize variables.
filename = 'M_code.txt';
delimiter = '';

%% Format string for each line of text:
%   column1: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Create output variable
SymbolicOutput = [dataArray{1:end-1}];

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% replace brackets and commas
for lineIndex = 1:length(SymbolicOutput)
    SymbolicOutput(lineIndex) = regexprep(SymbolicOutput(lineIndex), '\(', '[');
    SymbolicOutput(lineIndex) = regexprep(SymbolicOutput(lineIndex), '\)', ']');
    SymbolicOutput(lineIndex) = regexprep(SymbolicOutput(lineIndex), '\,', '][');
end

%% Convert declarations
for lineIndex = 1:length(SymbolicOutput)
    str = char(SymbolicOutput(lineIndex));
    if ~isempty(regexp(str,'zeros', 'once'))
        index1 = regexp(str,' = zeros[','once')-1;
        index2 = regexp(str,' = zeros[','end','once')+1;
        index3 = regexp(str,'\]\[','once')-1;
        index4 = index3 + 3;
        index5 = max(regexp(str,'\]'))-1;
        str1 = {'float '};
        str2 = str(1:index1);
        str3 = '[';
        str4 = str(index2:index3);
        str5 = '][';
        str6 = str(index4:index5);
        str7 = '];';
        if isempty(regexp(str,'\[1\]\;', 'once'))
            SymbolicOutput(lineIndex) = strcat(str1,str2,str3,str4,str5,str6,str7);
        else
            SymbolicOutput(lineIndex) = strcat(str1,str2,str3,str4,str7);
        end
    end
end

%% Convert indexing
strLeft = '\[';
strRight = '\]';
for arrayIndex = 1:24
    strIndex = int2str(arrayIndex);
    strRep = sprintf('[%d]',(arrayIndex-1));
    strPat = strcat(strLeft,strIndex,strRight);
    for lineIndex = 1:length(SymbolicOutput)
        str = char(SymbolicOutput(lineIndex));
        if isempty(regexp(str,'float','once'))
            SymbolicOutput(lineIndex) = {regexprep(str, strPat, strRep)};
        end
    end
end

%% replace divisions
for lineIndex = 1:length(SymbolicOutput)
    SymbolicOutput(lineIndex) = regexprep(SymbolicOutput(lineIndex), '/2', '*0.5');
    SymbolicOutput(lineIndex) = regexprep(SymbolicOutput(lineIndex), '/4', '*0.25');
end

%% Write to file

fid = fopen('C_code.txt','wt');
for lineIndex = 1:length(SymbolicOutput)
    fprintf(fid,char(SymbolicOutput(lineIndex)));
    fprintf(fid,'\n');
end
fclose(fid);
