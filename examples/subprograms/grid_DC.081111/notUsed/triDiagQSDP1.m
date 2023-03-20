%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [QMat0,cVect0,lbd,ubd,K0,L0] = triDiagQSDP1(nDim,randSeed,fileSaveSW); 

%
% Nonconvex tridiagonal QSDP
%

% Input
%   nDim
%   randSeed
%   fileSaveSW 
%       = 1 ---> To save the output in the file 'arrowQSDP2_k.mat', 
%                   where k = num2str(kBlock). 
%       = 0 ---> Not to save the output in any file --- default. 

%
% MMat(WMat) \succeq \O, lbd <= x <= ubd 
%
% Here WMat = (1,x) (1,x)^T
%
% This quadratic iequality is represented as 
%   QMat0 vec(WMat) \succeq 0, lbd <= x <= ubd
% 
% [L0.s*L0.s,K0.s*K0.s] = size(QMat0) 
%
% The objective function is represented as 
%   cVect0' * vec(WMat)
% to be minimized
%

if nargin == 0
    nDim = 3;
    randSeed = 2008; 
    fileSaveSW = 0;
elseif nargin == 1
    randSeed = 2008; 
    fileSaveSW = 0;
elseif nargin == 2
    fileSaveSW = 0;
end    

rand('seed',randSeed); 

K0.s = 1+nDim; % the size of WMat
L0.s = nDim; % the size of MMat

% lbd = [];
% ubd = [];

lbd = -ones(nDim,1);
ubd = ones(nDim,1); 

L0.f = 0;
L0.l = 0;

QMat0 = [];
for i=1:nDim
    for j= 1:nDim
        if i==j 
            % i== j|
            % This part can be modified as long as oneCoefMat is K0.s*K0.s
            % symmetric matirces ---> 
            oneCoefMat = sparse(K0.s,K0.s); 
            oneCoefMat(1,1) = rand(1,1); 
            oneCoefMat(i+1,i+1) = -1;
            % <---
            QMat0 = [QMat0;reshape(oneCoefMat,1,K0.s*K0.s)];
        elseif j == i+1
            % |i-j| == 1
            % This part can be modified as long as oneCoefMat is K0.s*K0.s
            % symmetric matirces ---> 
            oneCoefMat = sparse(K0.s,K0.s); 
            ai = rand(1,1) - 0.5;
            bi = rand(1,1) - 0.5;
            ci = rand(1,1) - 0.5; 
            oneCoefMat(1,1+i) = -ai/2; 
            oneCoefMat(1,1+j) = -bi/2; 
            oneCoefMat(1+i,1) = -ai/2; 
            oneCoefMat(1+j,1) = -bi/2; 
                % nonconvexity ---> 
                oneCoefMat(1+i,1+j) = ci;
                oneCoefMat(1+j,1+i) = ci;
                % <--- nonconvexity 
            % <---
            QMat0 = [QMat0;reshape(oneCoefMat,1,K0.s*K0.s)];
        elseif i == j+1
            % |i-j| == 1
            % This part can not be modified to keek the symmetry of the
            % matrix inequality ---> 
            p = (j-1)*nDim+i; 
            QMat0 = [QMat0; QMat0(p,:)];
            % <---
        else
            % |i-j| > 1
            % This part can not be modified to keek the tridaigonality of the
            % matrix inequality ---> 
            QMat0 = [QMat0; sparse(1,K0.s*K0.s)];
            % <--- 
        end
    end
end

cMat = sparse(K0.s,K0.s);
dVect = rand(1,nDim)/2;
cMat(1,2:K0.s) = dVect;
cMat(2:K0.s,1) = dVect';
cVect0 = reshape(cMat,K0.s*K0.s,1); 

debugSW = 0;
if debugSW == 1
    K0.s
    L0.s
%    mDim = size(QMat0,1); 
%     for i=1:mDim
%         i
%         full(reshape(QMat0(i,:),K0.s,K0.s))        
%     end 
%    full(reshape(cVect0,K0.s,K0.s))
    full(cVect0')
    full(QMat0)
end
debugSW = 0;

if fileSaveSW == 1
    fileName = strcat('triDiagQSDP1_',num2str(nDim),'.mat'); 
    save(fileName,'QMat0','cVect0','lbd','ubd','K0','L0'); 
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
