function [A,b,c,K,J] = arrowTriDLOP(sDim,kBlock,randSeed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: sDim, nDim, randSeed
%   QOP:
%   minimize    c^T vec ((1,x)^T (1,x))
%   subject to  Q vec ((1,x)^T (1,x)) \in J,  lbd \leq x \leq ubd,  x \in R^sDim
%
% Convex arrow type QSDP
%   c^T vec ((1,x)^T (1,x)) --- linear 
%   Q vec ((1,x)^T (1,x)) --- nonconvex, arrow type
%   Q : nDim \times nDim, where nDim = 2*kBlock+1 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output --- SDP relaxation of the QSDP 
%
%   <CoLO> Conic form Linear Optimization Problem
%       minimize   c^T x 
%       subject to A x - b \in J,  x \in K
%   Here 
%       x : a column vector variable. 
%       K.f --- the number of free variables, e.g., K.f = [], 0 or 10
%       K.l --- the number of LP variables,e.g., K.l = [], 0 or 12 
%       K.q --- the structure of SOCP variables, e.g., K.q = [], 3 or [3,5] 
%       K.s --- the structure of SDP variables, e.g., K.s = [], 4 or [2,4]
%       J.f --- the number of equality constraints, e.g., J.f = [], 0 or 6
%       J.l --- the number of LP inequality constraints, e.g., J.l = [], 0 or 7
%       J.q --- the structure of SOCP constraints, e.g., J.q = [], 8 or [2,3] 
%       J.s --- the structure of SDP constraints, e.g., J.s = [], 2 or [3,6]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [Q,c,sDim,J,lbd,ubd] = arrowTriNonConvex(sDim,kBlock,randSeed); 

[Qt,c, sDim, Js] = mexArrowTriDQOP(sDim,kBlock,randSeed);

J.s = Js; 
lbd = [];
ubd = [];

[A,b,c,J,K] = QSDPtoSDPrelaxation(Qt',c,sDim,J,lbd,ubd);

return