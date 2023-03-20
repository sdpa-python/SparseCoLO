function [QMat0,cVect0,J,K0] = addLbdUbd(QMat0,cVect0,lbd,ubd,J0,K0); 

rowSizeQMat0 = size(QMat0,1);
% colSizeQMat0 = size(QMat0,2);

if isfield(J0,'f') && J0.f > 0
    J.f = J0.f; 
    QMat0Free = QMat0(1:J0.f,:); 
    QMat0 = QMat0(J0.f+1:rowSizeQMat0,:); 
    rowSizeQMat0 = rowSizeQMat0 - J0.f;
else
    QMat0Free = [];
end

if isfield(J0,'l') && J0.l > 0
    J.l = J0.l;
else
    J.l = 0; 
end

if isfield(J0,'s') && length(J0.s) > 0
    J.s = J0.s;
% else
%     J.s = 0; 
end

QMatBound = [];
for i = 1:K0.s-1
    % 0 <= (x_i -lbd(i))(ubd(i) -x_i) = -x_i^2 +(lbd(i)+ubd(i) x_i -lbd(i)*ubd(i)
    oneCoefMat = sparse(K0.s,K0.s);
    oneCoefMat(1,1) = -lbd(i)*ubd(i); % - lbd(i)+ubd(i)
    oneCoefMat(1,1+i) = 0.5*(lbd(i)+ubd(i)); % + 0.5(lbd(i)+ubd(i)) x_i
    oneCoefMat(1+i,1) = 0.5*(lbd(i)+ubd(i)); % + 0.5(lbd(i)+ubd(i)) x_i
    oneCoefMat(1+i,1+i) = -1;  %-x_i^2
    QMatBound = [QMatBound; reshape(oneCoefMat,1,K0.s*K0.s)]; 
end
J.l = J.l + K0.s-1; 

QMat0 = [QMat0Free; QMatBound; QMat0]; 

return