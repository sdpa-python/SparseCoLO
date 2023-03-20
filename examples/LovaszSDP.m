function [A,b,c,K,J] = LovaszSDP(type,noOfVertices,nodeDegree,randSeed); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
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
% Output --- SDP relaxation of the problem of computing Lovasz number
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

[statusSW,costMatrix,maxDegree,minDegree,aveDegree] = generateNetwork4(type,noOfVertices,nodeDegree,randSeed);

nDim = noOfVertices;
PMat = speye(nDim,nDim);
PMat(1:nDim-1,2:nDim) = PMat(1:nDim-1,2:nDim) - speye(nDim-1,nDim-1);
c = sparse(nDim*nDim,1);
c(nDim*nDim,1) = -1;
A = reshape(PMat*PMat',1,nDim*nDim);
b = 1;
for i=1:nDim
    nzIdx = find(costMatrix(i,:));
    for j=nzIdx
        AMat = sparse(nDim,nDim);
        AMat(i,j) = 1;
        AMat(j,i) = 1;
        AMat = PMat*AMat*PMat';
        A = [A;reshape(AMat,1,nDim*nDim)];
        b = [b;0];
    end
end
K.s = nDim;

J.f = size(A,1);

return
