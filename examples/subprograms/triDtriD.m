%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q,c,sDim,J,lbd,ubd] = triDtriD(sDim,nDim,randSeed,fileSaveSW); 

% ======================================================================
%   QOP:
%   minimize    c^T vec ((1,x)^T (1,x))
%   subject to  Q vec ((1,x)^T (1,x)) \in J,  lbd \leq x \leq ubd,  x \in R^n
% 
% ======================================================================

%
% Convex tridiagonal QSDP
%   c^T vec ((1,x)^T (1,x)) --- linear 
%   Q vec ((1,x)^T (1,x)) --- tridiagonal, convex
% 

% Input
%   sDim
%   nDim
%   randSeed
%   fileSaveSW 
%       = 1 ---> To save the output in the file 'arrowQSDP2_k.mat', 
%                   where k = num2str(kBlock). 
%       = 0 ---> Not to save the output in any file --- default. 
%

if nargin == 0
    sDim = 3
    nDim = 4;
    randSeed = 2008; 
    fileSaveSW = 0;
elseif nargin == 1
    nDim = 4; 
    randSeed = 2008; 
    fileSaveSW = 0;
elseif nargin == 2
    randSeed = 2008; 
    fileSaveSW = 0;
elseif nargin == 3
    fileSaveSW = 0;    
end    

rand('seed',randSeed); 

K.s = 1+sDim; % the size of WMat
J.s = nDim; % the size of MMat

lbd = []; % ones(1,sDim);
ubd = []; % zeros(1,sDim);

Q = [];
for i=1:nDim
    for j= 1:nDim
        if i==j 
            % i== j|
            % This part can be modified as long as oneCoefMat is K.s*K.s
            % symmetric matirces ---> 
            oneCoefMat = sparse(K.s,K.s); 
            oneCoefMat(1,1) = 1;
            oneCoefMat(2:K.s,2:K.s) = -diag(rand(1,K.s-1)); 
            % <---
            Q = [Q;reshape(oneCoefMat,1,K.s*K.s)];
        elseif j == i+1
            % |i-j| == 1
            % This part can be modified as long as oneCoefMat is K.s*K.s
            % symmetric matirces ---> 
            oneCoefMat = sparse(K.s,K.s);
            oneCoefMat(1,2:K.s) = rand(1,K.s-1) - 0.5*ones(1,K.s-1); 
            oneCoefMat(2:K.s,1) = oneCoefMat(1,2:K.s)'; 
            % <---
            Q = [Q;reshape(oneCoefMat,1,K.s*K.s)];
        elseif i == j+1
            % |i-j| == 1
            % This part can not be modified to keek the symmetry of the
            % matrix inequality ---> 
            p = (j-1)*nDim+i; 
            Q = [Q; Q(p,:)];
            % <---
        else
            % |i-j| > 1
            % This part can not be modified to keek the tridaigonality of the
            % matrix inequality ---> 
            Q = [Q; sparse(1,K.s*K.s)];
            % <--- 
        end
    end
end
cMat = sparse(K.s,K.s);
dVect = (rand(1,sDim)-0.5*ones(1,sDim))/2;
cMat(1,2:K.s) = dVect;
cMat(2:K.s,1) = dVect';
c = reshape(cMat,K.s*K.s,1); 

debugSW = 0;
if debugSW == 1
    K.s
%   J.l
    J.s
%    mDim = size(Q,1); 
%     for i=1:mDim
%         i
%         full(reshape(Q(i,:),K.s,K.s))        
%     end 
%    full(reshape(c,K.s,K.s))
    full(c')
    full(Q)
end
debugSW = 0;

if fileSaveSW == 1
    fileName = strcat('triDtriD_',num2str(sDim),'_',num2str(nDim),'.mat'); 
    save(fileName,'Q','c','sDim','J','lbd','ubd'); 
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
