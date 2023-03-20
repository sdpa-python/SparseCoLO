function solveSDPs(outFileName,probIdxSet)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of SparseCoLO 
% Copyright (C) 2009 
% Masakazu Kojima Group
% Department of Mathematical and Computing Sciences
% Tokyo Institute of Technology
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

problem{1} = 'matFormat/arch8.mat';
problem{2} = 'matFormat/att532.mat';
problem{3} = 'matFormat/M1000.05.mat';
problem{4} = 'matFormat/M1000.10.mat';
problem{5} = 'matFormat/M1000.15.mat';
problem{6} = 'matFormat/maxG11.mat';
problem{7} = 'matFormat/maxG32.mat';
problem{8} = 'matFormat/mcp500-1.mat';
problem{9} = 'matFormat/qpG11.mat';
problem{10} = 'matFormat/thetaG11.mat';
problem{11} = 'matFormat/truss8.mat';
problem{12} = 'matFormat/sensor2Dim1000D.mat';
problem{13} = 'matFormat/sensor2Dim4000D.mat';
problem{14} = 'matFormat/sensor2Dim6000D.mat';
problem{15} = 'maxCutSDP(1,250,4,2009)'; 
problem{16} = 'maxCutSDP(1,500,4,2009)'; 
problem{17} = 'maxCutSDP(1,1000,4,2009)'; 
problem{18} = 'normMinSDP2(250,10,10,2009)'; 
problem{19} = 'normMinSDP2(500,10,10,2009)'; 
problem{20} = 'normMinSDP2(1000,10,10,2009)'; 
problem{21} = 'matFormat/mcp250-1.mat';
problem{22} = 'matFormat/ss30.mat';
problem{23} = 'matFormat/theta3.mat';
if nargin == 0
    outFileName = [];
    probIdxSet = [1:5,8:11,21:23];
elseif nargin == 1
    probIdxSet = [1:5,8:11,21:23];
end

%%%%%%%%%%
% SelectedParaSet = [parCoLO.domain,parCoLO.range,parCoLO.EQorLMI; parCoLO.domain,parCoLO.range,parCoLO.EQorLMI, ... 
%%%%%%%%%%
% SelectedParaSet = [0,0,1; 1,0,1; 2,0,2];
SelectedParaSet = [0,0,1; 1,0,1];
% SelectedParaSet = [0,0,1];
% SelectedParaSet = [0,0,2; 0,1,2];

for pp=probIdxSet; 
    dotDatsPosition0 = strfind(problem{pp},'.mat');
    if ~isempty(dotDatsPosition0)
        fileName = problem{pp}(1:dotDatsPosition0-1);
        fileName = problem{pp}(11:dotDatsPosition0-1);
        S = load(problem{pp}); 
        S.J.f = size(S.A,1);
        experimentCoLO3(S.A,S.b,S.c,S.K,S.J,SelectedParaSet,fileName,outFileName);
        clear S
    else
        fileName = problem{pp};
        [A,b,c,K,J] = eval(problem{pp});
        experimentCoLO3(A,b,c,K,J,SelectedParaSet,fileName,outFileName); 
        clear A b c K J
    end
end
return

function experimentCoLO3(A,b,c,K,J,SelectedParaSet,fileName,outFileName);

percent = char(37);
backS = char(92); % char(165);
aand = char(38);
perc = char(37);
lbrace = char(123);
rbrace = char(125);
lbracket = char(91);
rbracket = char(93);
pipe = char(124);

if nargin == 1
    fileName = A; 
    % Check whether ifileName has the extension '.mat' --->
    dotMatPosition0 = strfind(fileName,'.mat');
    % <--- Check whether fileName has the extension '.mat'
    if isempty(dotMatPosition0)
        [A,b,c,K,J] = eval(fileName);
    else
        S = load(fileName,'-mat');
        A = S.A;
        b = S.b;
        c = S.c;
        K = S.K;
        if isfield(S,'J')
            J = S.J;
        else
            J.f = size(A,1);
        end
        clear S
    end
% else
%     fileName = []; 
end

parCoLO.sedumipar.free = 1;
parCoLO.sedumipar.eps = 1.0e-9;
parCoLO.sedumipar.fid = 0;

parCoLO.sdpaOPTION.epsilonStar  = 1.0E-7;
parCoLO.sdpaOPTION.epsilonDash  = 1.0E-7;
parCoLO.sdpaOPTION.print  = '';

parCoLO.sdpt3OPTIONS.printlevel = 0;
% parCoLO.OPTIONSsdpt3.smallblkdim = 15;

% parCoLO.OPTIONsdpa.aggConeSize = 10; % [];

aggConeSizeFactor = 2; 

% parCoLO.SDPsolver = 'matFormat/sedumi'; 
% parCoLO.SDPsolver = 'matFormat/sdpa'; 

if ~isempty(outFileName)
    outFileID = fopen(outFileName,'a');
else
    outFileID = []; 
end

for solver = 1:3 
    if solver == 1
        parCoLO.SDPsolver = 'sedumi';
    elseif solver == 2
        parCoLO.SDPsolver = 'sdpa';
    else
        parCoLO.SDPsolver = 'sdpt3'; 
    end
    if solver == 1
        qq = 0;
    else
        qq = 0;
    end
    for q=0:qq
        if q == 0
            parCoLO.OPTIONsdpa.aggConeSize = [];
        else
            parCoLO.OPTIONsdpa.aggConeSize = (q+1)*aggConeSizeFactor;
        end
        
        SW = [];
        conversionTime = [];
        sdpCpuTime = [];
        primalObjectiveValues = [];
        pdGap = [];
        pfeasibility = [];
        dfeasibility = [];
        matAinfo = [];
        SDPinfo = [];
        coSpMatinfo = [];
        
        kk = size(SelectedParaSet,1);
        for k = 1:kk
            parCoLO.domain = SelectedParaSet(k,1);
            parCoLO.range = SelectedParaSet(k,2);
            parCoLO.EQorLMI = SelectedParaSet(k,3);
            %    fprintf('Start Experiment: domain =%2d, range =%2d, EQorLMI =%2d\n',...
            %            parCoLO.domain,parCoLO.range,parCoLO.EQorLMI);
            [x,y,infoCoLO,cliqueDomain,cliqueRange,LOP] = sparseCoLO(A,b,c,K,J,parCoLO);
            convTime = infoCoLO.CPUdomain + infoCoLO.CPUrange + infoCoLO.CPUEQorLMI;
            %
            [primalObjValue, dualObjValue, primalfeasibility, dualfeasibility] = evaluateCoLO(x,y,A,b,c,K,J,cliqueDomain,cliqueRange);
            SW = [SW; [parCoLO.domain, parCoLO.range, parCoLO.EQorLMI]];
            conversionTime = [conversionTime, convTime];
            sdpCpuTime = [sdpCpuTime, infoCoLO.CPUsolver];
            primalObjectiveValues = [primalObjectiveValues, primalObjValue];
            pfeasibility = [pfeasibility, primalfeasibility];
            dfeasibility = [dfeasibility, dualfeasibility];
            %
            if isfield(LOP.J,'f') && LOP.J.f == size(LOP.A,1)
                pdGap = [pdGap, abs(primalObjValue-dualObjValue)];
                [matrixA,SDPcone,coSpMat] = checkSparsityEQform(LOP.A,LOP.b,LOP.c,LOP.K);
            elseif isfield(LOP.K,'f') && LOP.K.f == size(LOP.A,2)
                pdGap = [pdGap, abs(dualObjValue-primalObjValue)];
                [matrixA,SDPcone,coSpMat] = checkSparsityEQform(-LOP.A',-LOP.c,-LOP.b,LOP.J);
            end
            %
            matAinfo = [matAinfo; [matrixA.size(1),matrixA.size(2) ,matrixA.nnz, matrixA.constSComp]];
            %
            SDPinfo = [SDPinfo; [SDPcone.noOfCones, SDPcone.sizeMax, SDPcone.sizeMin, SDPcone.sizeAve, SDPcone.volSDP]];
            %
            coSpMatinfo = [coSpMatinfo; [coSpMat.size, coSpMat.lbdNnz, coSpMat.nnz, coSpMat.nnzLMat]];
            fprintf('\n');
        end
        
        fprintf('\n');
        if ~isempty(fileName)
            fprintf('%s %s ',percent,fileName);
            if ~isempty(outFileID) 
                fprintf(outFileID,'%s %s ',percent,fileName);
            end
        end
        if strcmp(parCoLO.SDPsolver,'sedumi')
            fprintf('by SeDuMi\n');
            if ~isempty(outFileID)             
                fprintf(outFileID,'by SeDuMi\n');
            end
        elseif strcmp(parCoLO.SDPsolver,'sdpa')
            fprintf('by SDPAM with parCoLO.OPTIONsdpa.aggConeSize = %d\n',q*aggConeSizeFactor);
            if ~isempty(outFileID) 
                fprintf(outFileID,'by SDPAM with parCoLO.OPTIONsdpa.aggConeSize = %d\n',q*aggConeSizeFactor);
            end
        elseif strcmp(parCoLO.SDPsolver,'sdpt3')
            fprintf('by SDPT3 with parCoLO.OPTIONsdpa.aggConeSize = %d\n',q*aggConeSizeFactor);
            if ~isempty(outFileID) 
                fprintf(outFileID,'by SDPT3 with parCoLO.OPTIONsdpa.aggConeSize = %d\n',q*aggConeSizeFactor);
            end
        end
        fprintf('%sslover    parCoLO  |%scpu time |%s   matrix A          |%s Schur complement|   SDP blocks|\n',...
            percent,blanks(5),blanks(7),blanks(1));
        fprintf('%s         d  r EQ/LMI  cpuC     cpuS          sizeA          nnzA      nnzS     nnzL   noBl  maxBl\n',...
            percent);
        %
        if ~isempty(outFileID)
            fprintf(outFileID,'%sslover    parCoLO  |%scpu time |%s   matrix A          |%s Schur complement|   SDP blocks|\n',...
                percent,blanks(6),blanks(7),blanks(1));
            fprintf(outFileID,'%s        d  r EQ/LMI  cpuC     cpuS          sizeA          nnzA      nnzS     nnzL   noBl  maxBl\n',...
                percent);
        end
        for i=1:length(sdpCpuTime)
            if strcmp(parCoLO.SDPsolver,'sedumi')
                fprintf('sedumi %s ',aand);                
            elseif strcmp(parCoLO.SDPsolver,'sdpa')
                fprintf('sdpa   %s ',aand); 
            elseif strcmp(parCoLO.SDPsolver,'sdpt3')
                fprintf('sdpt3  %s ',aand); 
            end
            fprintf('%1d   %1d   %1d %s %6.1f %s %6.1f  %s %6d x%7d %s %7d %s %7d %s %6d %s %4d %s %4d%s%s\n',...
                SW(i,1),SW(i,2),SW(i,3),aand,conversionTime(i),aand,sdpCpuTime(i),...
                aand,matAinfo(i,1),matAinfo(i,2),aand,matAinfo(i,3),aand,coSpMatinfo(i,3),...
                aand,coSpMatinfo(i,4),aand,SDPinfo(i,1),aand,SDPinfo(i,2),backS,backS);
            %
            if ~isempty(outFileID)
                if strcmp(parCoLO.SDPsolver,'sedumi')
                    fprintf(outFileID,'sedumi %s',aand);
                elseif strcmp(parCoLO.SDPsolver,'sdpa')
                    fprintf(outFileID,'sdpa   %s',aand);
                elseif strcmp(parCoLO.SDPsolver,'sdpt3')
                    fprintf(outFileID,'sdpt3  %s',aand);
                end
                fprintf(outFileID,'%1d   %1d   %1d %s %6.1f %s %6.1f  %s %6d x%7d %s %7d %s %7d %s %6d %s %4d %s %4d%s%s\n',...
                    SW(i,1),SW(i,2),SW(i,3),aand,conversionTime(i),aand,sdpCpuTime(i),...
                    aand,matAinfo(i,1),matAinfo(i,2),aand,matAinfo(i,3),aand,coSpMatinfo(i,3),...
                    aand,coSpMatinfo(i,4),aand,SDPinfo(i,1),aand,SDPinfo(i,2),backS,backS);
            end
        end
        if length(sdpCpuTime) > 1
            fprintf('%s max primal objective value over the %2d cases = %+20.15e\n',percent,length(sdpCpuTime),max(primalObjectiveValues));
            fprintf('%s min primal objective value over the %2d cases = %+20.15e\n',percent,length(sdpCpuTime),min(primalObjectiveValues));
            fprintf('%s max primal obj. value- min primal obj. value = %+7.2e\n',percent,max(primalObjectiveValues)-min(primalObjectiveValues));
            %
            if ~isempty(outFileID)
                fprintf(outFileID,'%s max primal objective value over the %2d cases = %+20.15e\n',percent,length(sdpCpuTime),max(primalObjectiveValues));
                fprintf(outFileID,'%s min primal objective value over the %2d cases = %+20.15e\n',percent,length(sdpCpuTime),min(primalObjectiveValues));
                fprintf(outFileID,'%s max primal obj. value- min primal obj. value = %+7.2e\n',percent,max(primalObjectiveValues)-min(primalObjectiveValues));
            end
        end
        fprintf('\n');
        if ~isempty(outFileID) 
            fprintf(outFileID,'%s\n',percent);
        end
        % d = parCoLO.domain
        % r = parCoLO.range
        % EQ/LMI = parCoLO.EQorLMI
        % cpuC = the cpu time in second for conversion
        % cpuS = the cpu time in second for SeDuMi
        % gap = the gap between the primal and dual objective value
        % p.feas = the primal feasibility error
        % d.feas = the dual feasibility error
        % sizeA = the size of A
        % nnzA = the number of nonzeros in A
        % nnzS = the number of nonzeros in the Schur complement matrix
        % nnzL = the number of nonzeros in the sparse Cholesky factor of the
        %       Schur complement matrix
        % noBl = the number of SDP blocks
        % maxBl = the maximum size of SDP block
    end
end

return
    
