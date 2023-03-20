function [A,b,c,K,J] = fixLambdaLOP(sizePara,lambda,lbd,ubd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: sizePara,lambda,lbd,ubd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �����搶�C
% 
% ���x�����݂܂��񂪁C�⑫�ł��D
% 
% >���F
% >	\sum_{i} (A_i y_i + C_i y_i^{2}) 
% >		- \lambda \sum_{i} ( B_{i} y_ {i} +B_0) \succeq O 
% >�ɂ����āC
% >
% >�E�E�E�E�E�E�E
% 
% �Ə����܂������C�ړI�֐���
% 
% 	dll' * y	�̍ŏ���
% 
% �ł��D������l����ȊO�̐�������Ƃ��āCy_i �̔񕉐�����l���Ă��܂��D
% �܂Ƃ߂�ƁC
% 
% 	min.	\sum_{i} dll_i y_i
% 	s.t.	\sum_{i} (A_i y_i + C_i y_i^{2}) 
% 		- \lambda ( \sum_{i}  B_i y_ i +B_0) \succeq O 
% 		y_i \ge 0	(forall i)
% 
% �ł��D�i���̑O��block�Ίp���ł̎����̂悤�Ɂjy�̏�E��K���ɗ^��������
% ���C�H�w�I�ɂ͈Ӗ�������Ǝv���܂��D
% 
% --
% ���� �P��
% ������w��w�@ ��񗝍H�w�n������ �������w��U
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

