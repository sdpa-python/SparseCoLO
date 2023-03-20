function [A,b,c,K,J] = gPartitionSDP(type,noOfVertices,nodeDegree,randSeed); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% type
%   = 1 --- a sparse graph on [0,1] \times [0,1]
%   = 2 --- a sparse granp on the 2-dim unit sphere
%   = 3 --- a sparse graph on a torus
% noOfVertices
% nodeDegree >= 2
% randSeed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output --- SDP relaxation of a graph partitioning problem
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

if nodeDegree <= 1
    error('## nodeDegree needs to be not less than 2'); 
end

nodeDegree0 = nodeDegree;
nodeDegree = nodeDegree -1; 

[statusSW,costMatrix,maxDegree,minDegree,aveDegree] = ...
    generateNetwork4(type,noOfVertices,nodeDegree,randSeed);

nDim = noOfVertices;
PMat = speye(nDim,nDim);
PMat(1:nDim-1,2:nDim) = PMat(1:nDim-1,2:nDim) - speye(nDim-1,nDim-1);
costMatrix = (costMatrix + costMatrix');
d = costMatrix*ones(nDim,1);
costMatrix = costMatrix - diag(d);
costMatrix = PMat*costMatrix*PMat';
c = - reshape(costMatrix,nDim*nDim,1);
A = []; b = [];
for i=1:nDim
    AMat = sparse(nDim,nDim);
    AMat(i,i) = 1;
    AMat = PMat*AMat*PMat';
    A = [A;reshape(AMat,1,nDim*nDim)];
    b = [b; 0.25];
end
AMat = sparse(nDim,nDim);
AMat(nDim,nDim) = 1;
A = [A;reshape(AMat,1,nDim*nDim)];
b = [b; 0];
K.s = nDim;
J.f = size(A,1);

return
