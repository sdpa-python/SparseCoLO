function [A,b,c,K,J] = fixLambdaLOP(sizePara,lambda,lbd,ubd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: sizePara,lambda,lbd,ubd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    ubd = 1; 
end

[Q,c,sDim,J,lbd,ubd] = QSDPfixLambda(sizePara,lambda,lbd,ubd);

debugSW = 0; 
if debugSW == 1
    QOPinfo(Q,c,sDim,J);
    XXXXXXXXXX
end

[A,b,c,J,K] = QSDPtoSDPrelaxation(Q,c,sDim,J,lbd,ubd);

sDim

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function QOPinfo(Q,c,sDim,J);
fprintf('## sDim = %d, nDim = %d\n',sDim,J.s); 
Q = spones(Q); 
Qvect = sum(Q,1); 
spMat = reshape(Qvect,1+sDim,1+sDim);
figure(1);
spy(spMat);
Qvect = sum(Q,2); 
spMat = reshape(Qvect,J.s,J.s);
figure(2);
spy(spMat);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

