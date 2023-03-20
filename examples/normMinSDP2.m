
function [A,b,c,K,J] = normMinSDP2(q,r,t,randSeed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Norm minimization problem
%   minimize    \| F_{t+1} + \sum_{i=1}^r F_i z_i \|
% Here F_i : q \times r, q >> r.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% q : the number of rows of data matrices
% r : the number of columns of data matrices
% t : the number of data matrices
% randSeed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output --- a norm minimization problem
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

A = [];
b = [];
c = [];
K = [];

rand('seed',randSeed);

F_0 = 2*rand(q,r)-1.0;
for i=1:t
    F{i} = 2*rand(q,r)-1.0;
end
 
A = [];
for i=1:t
%    size([sparse(q,r),F{i}])
%    size(F{i}')
%    size([F{i}',sparse(r,q)])
        AMat = [[sparse(q,q),F{i}]; [ F{i}',sparse(r,r)]];       
        A = [A, reshape(AMat,(q+r)*(q+r),1)];
end

% AMat = speye(q+r,q+r);
% A = [reshape(AMat,(q+r)*(q+r),1),sparse((q+r)*(q+r),1),A];
% A = [sparse(1,2+t);A];
% A(1,2) = 1; 

% AMat = [[sparse(q,q),F_0]; [F_0',sparse(r,r)]];       
% b = [1; -reshape(AMat,(q+r)*(q+r),1)]; 
% c = sparse(2+t,1);
% c(1,1) = 1; 

% K.f = 1;
% K.q = 1+t;
% J.f = 1; 
% J.s = q+r; 

AMat = speye(q+r,q+r);
A = [reshape(AMat,(q+r)*(q+r),1),A];
AMat = [[sparse(q,q),F_0]; [F_0',sparse(r,r)]];       
b = [-reshape(AMat,(q+r)*(q+r),1)]; 
c = sparse(1+t,1);
c(1,1) = 1; 

K.f = t+1;
J.s = q+r; 

return