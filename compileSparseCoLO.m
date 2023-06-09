% complileSparseCoLO.m

%% Compiling Mex Files

% Check the existence of 'mex' subdirectory
currentDir = dir('*');
foundMexDir = 0;
currentDirNumber = length(currentDir);
for i=1:currentDirNumber
    if strcmp(currentDir(i).name,'mex')
        foundMexDir = 1;
        break;
    end
end
if foundMexDir == 0
    fprintf('Current directory is = %s\n',pwd);
    fprintf('Cannot find sub-directory "mex" \n');
    fprintf('Execute this script at the top directory of SparseCoLO\n');
    return;
end


eval('cd mex');

% 2013/03/15 Kojima <--- Hongbo Dong' comment 
% (1) Change lines 25-26 to:
% MLVerStr = regexp(version,'\d+','match');
% 
% MLVer = str2double(MLVerStr);
% ---> 
% MLVerStr = version;
% MLVer = str2num(MLVerStr(1:3));
MLVerStr = regexp(version,'\d+','match');
MLVer = str2double(MLVerStr);
% <--- 2013/03/15 Kojima <--- Hongbo Dong comment 

% if strcmp(computer, 'GLNXA64') && MLVer > 7.2
%     MexFlags = [' -O -largeArrayDims' ...
%                 ' CXXFLAGS="-fPIC -ansi -D_GNU_SOURCE" '];
% elseif strcmp(computer, 'GLNX86') && MLVer > 7.2
%     MexFlags = [' -O -largeArrayDims' ...
%                 ' CXXFLAGS="-fPIC -ansi -D_GNU_SOURCE -m32" '];
% else % Mac, Windows or Solaris
%      % MexFlags = ' -O -Dlinux=0 ';
%         MexFlags = ' -O ';
% end

% 2013/03/15 Kojima <--- Hongbo Dong's comment 
% 
% (2) Change all three "MLVer > 7.2" to
%
% (MLVer(1) > 7 || (MLVer(1) == 7&& MLVer(2) > 2))
% ---> 
% if strcmp(computer, 'GLNXA64') && MLVer > 7.2
if strcmp(computer, 'GLNXA64') && (MLVer(1) > 7 || (MLVer(1) == 7&& MLVer(2) > 2))
    MexFlags = [' -O -largeArrayDims' ...
                ' CXXFLAGS="-fPIC -ansi -D_GNU_SOURCE" '];
% elseif strcmp(computer, 'GLNX86') && MLVer > 7.2
elseif strcmp(computer, 'GLNX86') && (MLVer(1) > 7 || (MLVer(1) == 7&& MLVer(2) > 2))
    MexFlags = [' -O -largeArrayDims' ...
                ' CXXFLAGS="-fPIC -ansi -D_GNU_SOURCE -m32" '];
% elseif strcmp(computer, 'MACI64') && MLVer > 7.2
elseif strcmp(computer, 'MACI64') && (MLVer(1) > 7 || (MLVer(1) == 7&& MLVer(2) > 2))
    MexFlags = [' -O -largeArrayDims '];   
else % Mac, Windows or Solaris
     % MexFlags = ' -O -Dlinux=0 ';
        MexFlags = ' -O ';
end
% <---

% LIBfiles = ' conversion.cpp spvec.cpp polynomials.cpp sup.cpp ';
LIBfiles = [' ccputime.cpp'];
if ispc % Windows family create .obj files
        OBJfiles = strrep(LIBfiles,'.cpp','.obj');
else
        OBJfiles = strrep(LIBfiles,'.cpp','.o');
end

fprintf('Compiling Libraries...');
command = ['mex -c ' MexFlags LIBfiles];
eval(command);
fprintf('done\n');

mexFiles{1} = 'mexForestConvert.cpp';
mexFiles{2} = 'mexMaxSpanningTree2.cpp';
mexFiles{3} = 'mexPrimalOneSDP2.cpp';
mexFiles{4} = 'mexArrowTriDQOP.cpp';
mexFiles{5} = 'mexDiagTriDQOP.cpp';
for i=1:length(mexFiles)
    mexFileName = mexFiles{i};
    fprintf('Compiling %s...',mexFileName);
    command = ['mex ' MexFlags mexFileName OBJfiles];
    eval(command);
    fprintf('done\n');
end

eval('cd ../');
fprintf('Compilation finished successfully.\n');

