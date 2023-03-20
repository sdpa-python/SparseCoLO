function [A,b,c,J,K] = QSDPtoSDPrelaxation(Q,c,sDim,J,lbd,ubd);
%
%   QOP:
%   minimize    c^T vec ((1,x)^T (1,x))
%   subject to  Q vec ((1,x)^T (1,x)) \in J,  lbd \leq x \leq ubd,  x \in R^n
% 
% ======================================================================

if nargin == 6
    sparseSW.domain = 1; % 
    sparseSW.range = 1;
    sparseSW.EQorLMI = 1;
end

% Construction of an SDP relaxation problem ---> 
K.s = 1+sDim; 
if ~isempty(lbd) && ~isempty(ubd)
    % Incorporating lower and upper bounds into Q and J
    [A,c,J,K] = addLbdUbd(Q,c,lbd,ubd,J,K); 
elseif ~isempty(lbd)
    % Incorporating lower bounds into Q and J
    [A,c,J,K] = addLbd(Q,c,lbd,ubd,J,K); 
elseif ~isempty(ubd)
    % incorporating upper bounds into Q and J
    [A,c,J,K] = addUbd(Q,c,lbd,ubd,J,K); 
else
    % lbd = -\infty and ubd = +\infty 
    A = Q;
end
% fixing X_{00} = 1
[rowSizeA,colSizeA] = size(A); 
A = [speye(1,colSizeA); A]; 
b = speye(1+rowSizeA,1); 
if isfield(J,'f') && ~isempty(J.f)
    J.f = J.f + 1;
else
    J.f = 1;
end
% <--- Construction of an SDP relaxation problem

return


