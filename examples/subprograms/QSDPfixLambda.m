% function [QMat0,cVect0,lbd,ubd,K0,L0] = QSDPfixLambda(sizePara,lambda,lbd,ubd,fileSaveSW);
function [QMat0,cVect0,sDim,J,lbdVect,ubdVect] = QSDPfixLambda(sizePara,lambda,lbd,ubd);

fileSaveSW = 0; 

% 小島先生，
% 
% 何度もすみませんが，補足です．
% 
% >問題：
% >	\sum_{i} (A_i y_i + C_i y_i^{2}) 
% >		- \lambda \sum_{i} ( B_{i} y_ {i} +B_0) \succeq O 
% >において，
% >
% >・・・・・・・
% 
% と書きましたが，目的関数は
% 
% 	dll' * y	の最小化
% 
% です．半正定値制約以外の制約条件として，y_i の非負制約を考えています．
% まとめると，
% 
% 	min.	\sum_{i} dll_i y_i
% 	s.t.	\sum_{i} (A_i y_i + C_i y_i^{2}) 
% 		- \lambda ( \sum_{i}  B_i y_ i +B_0) \succeq O 
% 		y_i \ge 0	(forall i)
% 
% です．（この前のblock対角化での実験のように）yの上界を適当に与えた問題で
% も，工学的には意味があると思います．
% 
% --
% 寒野 善博
% 東京大学大学院 情報理工学系研究科 数理情報学専攻
% <kanno[at]mist.i.u-tokyo.ac.jp>

    
fileName = strcat('grid_DC',num2str(sizePara),'.mat'); 

S = load(fileName,'-mat','dll','AMat','BMat','CMat','BMat0'); 

nm = size(S.AMat,2);
% nm
L0.s = size(S.AMat{1},1);
%%%%%%%%%%%
J.s = L0.s; 
%%%%%%%%%%
% L0.s
% XXXXX
K0.s = 1+nm; % variables are y_1,... 
%%%%%%%%%% 
sDim = nm;
%%%%%%%%%%

bdConst = 1; 
if nargin == 0
    sizePara = 2; 
    lambda = 1000;
    lbd = 0; 
    ubd = 1; 
elseif nargin == 1
    lambda = 1000;
    lbd = 0; 
    ubd = 1; 
elseif nargin == 2
    lbd = 0; 
    ubd = 1; 
elseif nargin == 3
    ubd = 0; 
end

lbdVect = lbd*ones(1,nm);
ubdVect = ubd*ones(1,nm);

if isempty(lambda)
    lambda = 1000;
end
if isempty(lbd)
    lbd = zeros(1,nm); 
end
if isempty(ubd)
    ubd = ones(1,nm); 
end
if isempty(fileSaveSW)
    fileSaveSW = 0;
end
    
activeIdx = [1:nm];

% scaling ---> 
% lambda = 1.0e4 mu, lambda = 15,000 ---> mu = 1.5; 
% y_i = 1.0e-2 z_i, y_i = 0.1 ---> z_i = 10; 
for i=1:nm
    S.AMat{i} = 1.0e-2*S.AMat{i};
end
for i=1:nm
     S.BMat{i} = 1.0e-2*S.BMat{i};
end
for i=1:nm
    S.CMat{i} = 1.0e-4*S.CMat{i};
end
% S.BMat0 = 1.0e4*S.BMat0;
% <---

%%%%%%%
% mu = 1; % lambda = 15,000
%%%%%%%

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
objMat = sparse(K0.s,K0.s); 
objMat(2:K0.s,1) = S.dll(activeIdx,:)/2;
objMat(1,2:K0.s) = S.dll(activeIdx,:)'/2;
% objMat(1,K0.s) = -1;
% objMat(K0.s,1) = -1;
% full(objMat); 
cVect0 = reshape(objMat,K0.s*K0.s,1); 
clear objMat

% modifiedSW = 1;
% if modifiedSW == 0
%     QMat0 = [];
%     for p=1:L0.s
%         for q=1:L0.s
%             if p <= q
%                 oneCoefMat = sparse(K0.s,K0.s);
%                 for i=1:nm
%                     oneCoefMat(1,1+i) = S.AMat{i}(p,q) - lambda * S.BMat{i}(p,q);
%                     oneCoefMat(1+i,1) = oneCoefMat(1,1+i);
%                     oneCoefMat(1+i,1+i) = S.CMat{i}(p,q);
%                 end
%                 oneCoefMat(1,1) = - lambda * S.BMat0(p,q);
%                 QMat0 = [QMat0; reshape(oneCoefMat,1,K0.s*K0.s)];
%             else
%                 r = (q-1)*L0.s + p;
%                 QMat0 = [QMat0; QMat0(r,:)];
%             end
%         end
%     end
% else
    QMat0 = sparse(L0.s*L0.s,K0.s*K0.s);
    colPointer1 = 1; 
    QMat0(:,colPointer1) = -lambda*reshape(S.BMat0,L0.s*L0.s,1); 
    for i=1:nm
        colPointer1 = 1+i; 
        QMat0(:,colPointer1) = reshape(S.AMat{i},L0.s*L0.s,1) - lambda*reshape(S.BMat{i},L0.s*L0.s,1); 
        colPointer2 = i*K0.s + 1; 
        QMat0(:,colPointer2) = QMat0(:,colPointer1); 
        cloPointer1 = i*K0.s + 1+i; 
        QMat0(:,cloPointer1) = reshape(S.CMat{i},L0.s*L0.s,1);         
    end
% end

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
    fileName = strcat('QSDPfixLambda_',num2str(sizePara),'.mat'); 
    save(fileName,'QMat0','cVect0','lbd','ubd','K0','L0'); 
end

return

