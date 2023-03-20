%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [QMat0,cVect0,lbd,ubd,K0,L0] = arrowQSDP2(kBlock,randSeed,fileSaveSW); 

%
% Nonconvex arrow type QSDP
% 

% Input
%   kBlock
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
% The objective function is represented as 
%   cVect0' * vec(WMat)
% to be minimized
%
% K0.s = 1+2*kBlock+1; 
% L0.s = 2*kBlock+1; 
% 

if nargin == 0
    kBlock = 3;
    randSeed = 2008; 
    fileSaveSW = 0;
elseif nargin == 1
    randSeed = 2008; 
    fileSaveSW = 0;
elseif nargin == 2
    fileSaveSW = 0;
end    
    
rand('seed',randSeed); 

K0.f = 0; 

nDim = 2*kBlock+1; 
K0.s = 1+2*kBlock+1; % the size of WMat
L0.s = 2*kBlock+1; 

%
% lbd = [];
% ubd = [];
%
lbd = -ones(nDim,1);
ubd = ones(nDim,1); 
% ubd = zeros(nDim,1);

% the size of MMat
% [L0.s*L0.s,K0.s*K0.s] = size(QMat0) 
% 

% QMat0LP = [];

% for i=1:nDim
%     oneCoefMat = sparse(K0.s,K0.s);
%     oneCoefMat(1,1+i) = -0.5;
%     oneCoefMat(1+i,1) = -0.5;
%     QMat0LP = [QMat0LP;reshape(oneCoefMat,1,K0.s*K0.s)];
% end
% L0.l = nDim; 

QMat0 = [];
for i=1:nDim-1
    if mod(i,2) == 1
        for j=1:nDim-1
            if i==j
                % i== j|
                % This part can be modified as long as oneCoefMat is K0.s*K0.s
                % symmetric matirces --->
                % rand(1,1) - x_i^2
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
                oneCoefMat(1,1+i) = -ai/2;
                oneCoefMat(1,1+j) = -bi/2;
                oneCoefMat(1+i,1) = -ai/2;
                oneCoefMat(1+j,1) = -bi/2;
                % <---
                QMat0 = [QMat0;reshape(oneCoefMat,1,K0.s*K0.s)];
            else
                % |i-j| > 1
                % This part can not be modified to keek the tridaigonality of the
                % matrix inequality --->
                QMat0 = [QMat0; sparse(1,K0.s*K0.s)];
                % <---
            end
        end
        j = nDim;
        oneCoefMat = sparse(K0.s,K0.s);
        ci = rand(1,1) - 0.5; 
        oneCoefMat(1+i,1+j) = ci;
        oneCoefMat(1+j,1+i) = ci;
        QMat0 = [QMat0;reshape(oneCoefMat,1,K0.s*K0.s)];
    else
        for j=1:nDim-1
            if i==j
                % i== j|
                % This part can be modified as long as oneCoefMat is K0.s*K0.s
                % symmetric matirces --->
                oneCoefMat = sparse(K0.s,K0.s);
                oneCoefMat(1,1) = rand(1,1);
                oneCoefMat(i+1,i+1) = -1;
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
        % nonconvex term
        j = nDim;
        oneCoefMat = sparse(K0.s,K0.s);
        ci = rand(1,1) - 0.5; 
        oneCoefMat(1+i,1+j) = ci; % ci x_i*x_j
        oneCoefMat(1+j,1+i) = ci; % ci x_i*x_j 
        QMat0 = [QMat0;reshape(oneCoefMat,1,K0.s*K0.s)];
    end
end

%
i = nDim;
for j=1:nDim-1
    p = (j-1)*nDim+i;
    QMat0 = [QMat0; QMat0(p,:)];
end
%
i = nDim; j = nDim; 
% i== j|
% This part can be modified as long as oneCoefMat is K0.s*K0.s
% symmetric matirces --->
oneCoefMat = sparse(K0.s,K0.s);
oneCoefMat(1,1) = rand(1,1);
oneCoefMat(i+1,i+1) = -1;
QMat0 = [QMat0;reshape(oneCoefMat,1,K0.s*K0.s)];
%

cMat = sparse(K0.s,K0.s);
dVect = rand(1,nDim)/2;
cMat(1,2:K0.s) = dVect;
cMat(2:K0.s,1) = dVect';
cVect0 = reshape(cMat,K0.s*K0.s,1); 

% QMat0 = [QMat0LP;QMat0];

debugSW = 0;
if debugSW == 1
    K0.s
    L0.s
    for i=1:L0.s*L0.s
        i
        full(reshape(QMat0(i,:),K0.s,K0.s))        
    end
    full(reshape(cVect0,K0.s,K0.s))
end

debugSW = 0; 
if debugSW == 1
    fprintf('## Checking whether QMat0 is symmetric\n');
    errorSym = 0.0; 
    for i=1:L0.s-1
        for j=i+1:L0.s
            p = (i-1)*L0.s+j; 
            q = (j-1)*L0.s+i; 
            a = norm(QMat0(p,:) - QMat0(q,:),inf); 
            errorSym = errorSym + a; 
            if a > 1.0e-8
                fprintf('norm(QMat0(%d,:) - QMat0(%d,:),inf) > 1.0e-8\n',p,q);
            end
        end
    end
    fprintf('   Total error = %7.2e\n',errorSym);
    fprintf('## Checking whether each QMat0(i,:) is symmetric\n');
    for i=1:L0.s
        for j=i:L0.s
            p = (i-1)*L0.s+j; 
            tempMat = reshape(QMat0(p,:),K0.s,K0.s);
            a = norm(tempMat - tempMat',inf); 
            errorSym = errorSym + a; 
            if a > 1.0e-8
                fprintf('norm(QMat0(%d,:) - QMat0(%d,:)^T,inf) > 1.0e-8\n',p,p);
            end
        end
    end
    fprintf('   Total error = %7.2e\n',errorSym);    
end 

if fileSaveSW == 1
    fileName = strcat('arrowQSDP2_',num2str(kBlock),'.mat'); 
    save(fileName,'QMat0','cVect0','lbd','ubd','K0','L0'); 
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
