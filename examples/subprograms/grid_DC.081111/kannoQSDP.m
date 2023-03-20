function [QMat,cVect,K,L] = kannoQSDP(sizePara);

fileName = strcat('grid_DC',num2str(sizePara),'.mat'); 

S = load(fileName,'-mat','dll','AMat','BMat','CMat','BMat0'); 

nm = size(S.AMat,2);
L.s = size(S.AMat{1},1);
K.s = 1+size(S.dll,1)+1; % variables are y_1,... and \lambda

% S.dll = S.dll; % + 1.0e-3*rand(nm,1);

% scaling ---> 
% \lambda = 1.0e4 \mu 
% y_i = 1.0e-2 z_i 
for i=1:nm
    S.AMat{i} = 1.0e-2*S.AMat{i};
end
for i=1:nm
    S.BMat{i} = 1.0e2*S.BMat{i};
end
for i=1:nm
    S.CMat{i} = 1.0e-4*S.CMat{i};
end
S.BMat0 = 1.0e4*S.BMat0;
% <---

debugSW = 0;
if debugSW == 1
    for i=1:nm
        fprintf('norm(S.AMat{%2d},inf) = %7.2e\n',i,norm(S.AMat{i},inf));
    end
    for i=1:nm
        fprintf('norm(S.BMat{%2d},inf) = %7.2e\n',i,norm(S.BMat{i},inf));
    end
    for i=1:nm
        fprintf('norm(S.CMat{%2d},inf) = %7.2e\n',i,norm(S.CMat{i},inf));
    end
    fprintf('norm(S.BMat%d,inf) = %7.2e\n',0,norm(S.BMat0,inf));
end

% 
objMat = sparse(K.s,K.s); 
objMat(2:K.s-1,1) = S.dll/2;
objMat(1,2:K.s-1) = S.dll'/2;
% full(objMat); 
cVect = reshape(objMat,K.s*K.s,1); 
clear objMat

QMatLP = [];
% Nonnegativity condition 
for i = 1:nm+1
    oneCoefMat = sparse(K.s,K.s);
    oneCoefMat(1,1+i) = 0.5;
    oneCoefMat(1+i,1) = 0.5;
     if i == nm+1
         oneCoefMat(1,1) = -0.72; 
     end
    QMatLP = [QMatLP; reshape(oneCoefMat,1,K.s*K.s)]; 
end
aVect = 2*ones(1,nm+1); 
for i = 1:nm+1
    oneCoefMat = sparse(K.s,K.s);
    oneCoefMat(1,1+i) = -0.5;
    oneCoefMat(1+i,1) = -0.5;
    oneCoefMat(1,1) = aVect(i); 
    QMatLP = [QMatLP; reshape(oneCoefMat,1,K.s*K.s)]; 
end
aVect = 2*ones(1,nm+1); 
for i = 1:nm+1
    oneCoefMat = sparse(K.s,K.s);
    oneCoefMat(1+i,1+i) = -1;
    oneCoefMat(1,1+i) = aVect(i)/2; 
    oneCoefMat(1+i,1) = aVect(i)/2;     
    QMatLP = [QMatLP; reshape(oneCoefMat,1,K.s*K.s)]; 
end

aVect = 2*ones(1,nm+1); 
for i = 1:nm
    oneCoefMat = sparse(K.s,K.s);
    oneCoefMat(1+i,K.s) = -0.5;
    oneCoefMat(K.s,1+i) = -0.5;
    oneCoefMat(1,1+i) = aVect(i)/2; 
    oneCoefMat(1+i,1) = aVect(i)/2;     
    QMatLP = [QMatLP; reshape(oneCoefMat,1,K.s*K.s)]; 
end

aVect = 2*ones(1,nm+1); 
for i = 1:nm
    oneCoefMat = sparse(K.s,K.s);
    oneCoefMat(1+i,K.s) = -0.5;
    oneCoefMat(K.s,1+i) = -0.5;
    oneCoefMat(1,K.s) = aVect(i)/2; 
    oneCoefMat(K.s,1) = aVect(i)/2;     
    QMatLP = [QMatLP; reshape(oneCoefMat,1,K.s*K.s)]; 
end

L.l = 5*nm+3; 

QMat = [];
for p=1:L.s
    for q=1:L.s
        if p <= q
            oneCoefMat = sparse(K.s,K.s);
            for i=1:nm
                oneCoefMat(1,1+i) = S.AMat{i}(p,q);
                oneCoefMat(1+i,1) = oneCoefMat(1,1+i);
                oneCoefMat(1+i,1+i) = S.CMat{i}(p,q);
                oneCoefMat(K.s,1+i) = -S.BMat{i}(p,q);
                oneCoefMat(1+i,K.s) = oneCoefMat(K.s,1+i);
            end
            oneCoefMat(1,K.s) = -S.BMat0(p,q);
            oneCoefMat(K.s,1) = oneCoefMat(1,K.s); 
            QMat = [QMat; reshape(oneCoefMat,1,K.s*K.s)]; 
        else
            r = (q-1)*L.s + p; 
            QMat = [QMat; QMat(r,:)]; 
        end
    end
end

QMat = [QMatLP;QMat];

debugSW = 0;
if debugSW == 1
    K.s
    L.s
    for i=1:L.s*L.s
        i
        full(reshape(QMat(i,:),K.s,K.s))        
    end
    full(reshape(cVect,K.s,K.s))
end

debugSW = 0; 
if debugSW == 1
    fprintf('## Checking whether QMat is symmetric\n');
    errorSym = 0.0; 
    for i=1:L.s-1
        for j=i+1:L.s
            p = (i-1)*L.s+j; 
            q = (j-1)*L.s+i; 
            a = norm(QMat(p,:) - QMat(q,:),inf); 
            errorSym = errorSym + a; 
            if a > 1.0e-8
                fprintf('norm(QMat(%d,:) - QMat(%d,:),inf) > 1.0e-8\n',p,q);
            end
        end
    end
    fprintf('   Total error = %7.2e\n',errorSym);
    fprintf('## Checking whether each QMat(i,:) is symmetric\n');
    for i=1:L.s
        for j=i:L.s
            p = (i-1)*L.s+j; 
            tempMat = reshape(QMat(p,:),K.s,K.s);
            a = norm(tempMat - tempMat',inf); 
            errorSym = errorSym + a; 
            if a > 1.0e-8
                fprintf('norm(QMat(%d,:) - QMat(%d,:)^T,inf) > 1.0e-8\n',p,p);
            end
        end
    end
    fprintf('   Total error = %7.2e\n',errorSym);    
end            

return

