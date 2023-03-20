function [xVect,WMat] = solveSDPrelaxOfQSDP(Q,c,sDim,J,lbd,ubd,parCoLO);
%
%   QOP:
%   minimize    c^T vec((1,x)^T (1,x))
%   subject to  Q vec((1,x)^T (1,x)) \in J,  lbd \leq x \leq ubd,  x \in R^sDim
% 
%   Here 
%       x : a column vector variable with dimension sDim
%       J.f --- the number of equality constraints, e.g., J.f = [], 0 or 6
%       J.l --- the number of LP inequality constraints, e.g., J.l = [], 0 or 7
%       J.q --- the structure of SOCP constraints, e.g., J.q = [], 8 or [2,3] 
%       J.s --- the structure of SDP constraints, e.g., J.s = [], 2 or [3,6]
% ======================================================================
%
%   <parameters>
%
%   parCoLO.method 
%       = 0 --- do not apply any method
%       = 2 --- sedumi
%
%   parCoLO.parSeDuMi --- parameters for sedumi
%
%   parCoLO.domain
%       = 0 --- exploiting no sparsity in the domain space
%       = 1 --- applying dConvCliqueTree (default)
%       = 2 --- applying dConvBasisRep
%   This switch is irrevant if no psd matrix variable is involved.  
%
%   parCoLO.range
%       = 0 --- exploiting no sparsity in the range space 
%       = 1 --- applying rConvCliqueTree 
%       = 2 --- applying rConvMatDecomp (default) 
%   This switch is irrevant if no matrix inequality constraint is involved. 
%
%   parCoLO.EQorLMI
%       = 1 --- applying CoLOtoEQform to obtain an equality standard from, 
%               which can be solved by many existing software packages.
%               (default)
%       = 2 --- applying CoLOtoLMIform to obtain an LMI standard form, 
%               which can be solved by many existing software packages. 
%
%   One recommended choice of parameters: 
%       parCoLO.domain = 1;  % dConvCliqueTree  ---> equalities 
%       parCoLO.range = 2;   % rConvMatDecomp   ---> equalities 
%       parCoLO.EQorLMI = 1; % CoLOtoEQform     ---> equality standard form
%   The other recommended choice of parameters: 
%       parCoLO.domain = 2;  % dConvBasisRep    ---> matrix inequalities 
%       parCoLO.range = 1;   % rConvCliqueTree  ---> matrix inequalities 
%       parCoLO.EQorLMI = 2; % CoLOtoLMIform    ---> LMI standard form
% 
% ======================================================================

% < Sample excecution 1>
% >> parCoLO.domain = 2;
% >> parCoLO.range = 1;
% >> parCoLO.EQorLMI = 2;
% >> [Q,c,sDim,J,lbd,ubd] = triDtriD(6,6,2009);
% >> [xVect,clique,,WWMat] = solveSDPrelaxOfQSDP(Q,c,sDim,J,lbd,ubd,parCoLO);
% Applying the d-space conversion method using basis representation
% Applying the r-space conversion method using clique trees
% LOP to be converted into LMI standard form is already LMI standard form
% SeDuMi 1.1R3 by AdvOL, 2006 and Jos F. Sturm, 1998-2003.
% Alg = 2: xz-corrector, Adaptive Step-Differentiation, theta = 0.250, beta = 0.500
% Put 1 free variables in a quadratic cone
% eqs m = 13, order n = 21, dim = 63, blocks = 9
% nnz(A) = 91 + 0, nnz(ADA) = 169, nnz(L) = 91
%  it :     b*y       gap    delta  rate   t/tP*  t/tD*   feas cg cg  prec
%   0 :            1.51E+01 0.000
%   1 :   1.79E-01 5.09E+00 0.000 0.3375 0.9000 0.9000   1.40  1  1  7.0E+00
%   2 :   4.07E-01 1.37E+00 0.000 0.2698 0.9000 0.9000   2.63  1  1  1.1E+00
%   3 :   5.93E-01 1.06E-01 0.000 0.0770 0.9900 0.9900   1.14  1  1  7.6E-02
%   4 :   6.03E-01 2.23E-02 0.000 0.2113 0.9000 0.9000   1.02  1  1  1.6E-02
%   5 :   6.06E-01 6.33E-04 0.000 0.0283 0.9900 0.9900   1.00  1  1  4.5E-04
%   6 :   6.06E-01 7.50E-05 0.000 0.1186 0.9091 0.9000   1.00  1  1  5.8E-05
%   7 :   6.06E-01 9.81E-06 0.000 0.1307 0.9089 0.9000   1.00  1  1  8.0E-06
%   8 :   6.06E-01 4.01E-07 0.000 0.0409 0.9901 0.9900   1.00  1  1  3.5E-07
%   9 :   6.06E-01 6.05E-08 0.000 0.1508 0.9092 0.9000   1.00  1  1  5.9E-08
%  10 :   6.06E-01 3.80E-09 0.202 0.0628 0.9900 0.9902   1.00  1  1  3.5E-09
%  11 :   6.06E-01 3.99E-10 0.000 0.1050 0.9184 0.9000   1.00  2  2  5.0E-10
% 
% iter seconds digits       c*x               b*y
%  11      0.2  10.0  6.0625766579e-01  6.0625766573e-01
% |Ax-b| =   1.3e-10, [Ay-c]_+ =   1.3E-10, |x|=  1.4e+00, |y|=  1.6e+00
% 
% Detailed timing (sec)
%    Pre          IPM          Post
% 1.000E-02    2.500E-01    5.000E-02    
% Max-norms: ||b||=3.838117e-01, ||c|| = 1,
% Cholesky |add|=0, |skip| = 0, ||L.L|| = 6.05702.
% >> full(xVect')
% 
% ans =
% 
%     0.6982   -0.1336   -0.3701   -0.3313   -0.0652   -0.6158
% >> clique
% 
% clique = 
% 
%        NoC: 6
%     NoElem: [2 2 2 2 2 2]
%       maxC: 2
%       minC: 2
%        Set: {[1 2]  [1 3]  [1 4]  [1 5]  [1 6]  [1 7]}
% >> full(incidenceMatrixT) --- incidence matrix for a clique tree
% 
% ans =
% 
%      1     0     0     1     0
%     -1     0     0     0     0
%      0     0     1     0     1
%      0     0    -1     0     0
%      0     1     0    -1     0
%      0    -1     0     0    -1
% >> full(WWMat)
% 
% ans =
% 
%     1.0000    0.6982   -0.1336   -0.3701   -0.3313   -0.0652   -0.6158
%     0.6982    0.4875         0         0         0         0         0
%    -0.1336         0    0.0178         0         0         0         0
%    -0.3701         0         0    0.1370         0         0         0
%    -0.3313         0         0         0    0.1097         0         0
%    -0.0652         0         0         0         0    0.0043         0
%    -0.6158         0         0         0         0         0    0.3792

% <Sample excecution 2> 
% parCoLO.domain = 2;
% parCoLO.range = 1;
% parCoLO.EQorLMI = 2;
% >> [Q,c,sDim,J,lbd,ubd] = arrowTriNonConvex(20,20,2009);
% >> [xVect,clique,incidenceMatrixT,WWMat] = solveSDPrelaxOfQSDP(Q,c,sDim,J,lbd,ubd);
% Applying the d-space conversion method using clique trees
% Applying the r-space conversion method using matrix decomposition
% Conversion into an equality standard form
% SeDuMi 1.1R3 by AdvOL, 2006 and Jos F. Sturm, 1998-2003.
% Alg = 2: xz-corrector, Adaptive Step-Differentiation, theta = 0.250, beta = 0.500
% eqs m = 125, order n = 104, dim = 490, blocks = 23
% nnz(A) = 4569 + 0, nnz(ADA) = 15391, nnz(L) = 7758
%  it :     b*y       gap    delta  rate   t/tP*  t/tD*   feas cg cg  prec
%   0 :            5.39E+00 0.000
%   1 :  -3.14E-02 1.36E+00 0.000 0.2519 0.9000 0.9000   2.34  1  1  6.7E+00
%   2 :  -5.08E-01 4.04E-01 0.000 0.2976 0.9000 0.9000   1.53  1  1  1.5E+00
%   3 :  -4.33E-01 1.42E-01 0.000 0.3505 0.9000 0.9000   1.54  1  1  4.3E-01
%   4 :  -3.86E-01 4.87E-02 0.000 0.3435 0.9000 0.9000   1.08  1  1  1.5E-01
%   5 :  -3.62E-01 1.80E-02 0.000 0.3702 0.9000 0.9000   0.89  1  1  6.0E-02
%   6 :  -3.41E-01 8.11E-03 0.000 0.4502 0.9000 0.9000   0.71  1  1  3.2E-02
%   7 :  -3.25E-01 4.36E-03 0.000 0.5379 0.9000 0.9000   0.61  1  1  2.0E-02
%   8 :  -3.09E-01 2.29E-03 0.000 0.5256 0.9000 0.9000   0.49  1  1  1.3E-02
%   9 :  -2.95E-01 1.29E-03 0.000 0.5618 0.9000 0.9000   0.47  1  1  8.7E-03
%  10 :  -2.83E-01 7.06E-04 0.000 0.5476 0.9000 0.9000   0.53  1  1  5.5E-03
%  11 :  -2.75E-01 4.51E-04 0.000 0.6395 0.9000 0.9000   0.47  1  1  4.2E-03
%  12 :  -2.65E-01 2.25E-04 0.000 0.4976 0.9000 0.9000   0.58  1  1  2.4E-03
%  13 :  -2.58E-01 1.25E-04 0.000 0.5551 0.9000 0.9000   0.53  1  1  1.6E-03
%  14 :  -2.50E-01 4.96E-05 0.000 0.3977 0.9000 0.9000   0.71  1  1  7.0E-04
%  15 :  -2.47E-01 2.12E-05 0.000 0.4277 0.9000 0.9000   0.75  1  1  3.3E-04
%  16 :  -2.44E-01 6.50E-06 0.000 0.3065 0.9000 0.9000   0.88  1  1  1.1E-04
%  17 :  -2.44E-01 2.34E-06 0.000 0.3596 0.9000 0.9000   0.93  1  1  4.0E-05
%  18 :  -2.43E-01 6.23E-07 0.000 0.2665 0.9000 0.9000   0.97  1  1  1.1E-05
%  19 :  -2.43E-01 1.19E-07 0.000 0.1906 0.8679 0.9000   0.99  1  1  2.4E-06
%  20 :  -2.43E-01 9.58E-09 0.081 0.0806 0.9900 0.9900   1.00  1  1  2.0E-07
%  21 :  -2.43E-01 1.81E-09 0.000 0.1889 0.9000 0.9000   1.00  1  1  3.7E-08
%  22 :  -2.43E-01 8.32E-11 0.469 0.0460 0.9900 0.9900   1.00  2  2  1.7E-09
%  23 :  -2.43E-01 1.61E-11 0.000 0.1939 0.9000 0.9000   1.00  3  8  3.3E-10
% 
% iter seconds digits       c*x               b*y
%  23      0.7   Inf -2.4321681223e-01 -2.4321681216e-01
% |Ax-b| =   1.2e-09, [Ay-c]_+ =   5.2E-12, |x|=  5.9e+00, |y|=  4.8e+01
% 
% Detailed timing (sec)
%    Pre          IPM          Post
% 4.000E-02    7.000E-01    0.000E+00    
% Max-norms: ||b||=1, ||c|| = 9.761890e-01,
% Cholesky |add|=1, |skip| = 1, ||L.L|| = 500000.

if nargin <= 6
    parCoLO.method = 2; 
    parCoLO.domain = 1;  
    parCoLO.range = 2;
    parCoLO.EQorLMI = 1;
    parCoLO.parSeDuMi.fid = 1;
    parCoLO.parSeDuMi.free = 1;
    parCoLO.parSeDuMi.eps = 1.0e-9; 
end
if ~isfield(parCoLO,'method') 
    parCoLO.method = 2; 
end
if ~isfield(parCoLO,'parSeDuMi') 
    parCoLO.parSeDuMi.fid = 1;
    parCoLO.parSeDuMi.free = 1;
    parCoLO.parSeDuMi.eps = 1.0e-9; 
end

if nargin <= 4
    lbd = [];
    ubd = [];
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

% [x,y,sedumiInfo,cliqueDomain,cliqueRange] = sparseCoLO(A,b,c,K,J,conversionSW,pars); 
[x,y,infoCoLO,cliqueDomain,cliqueRange,LOP] = sparseCoLO(A,b,c,K,J,parCoLO); 

[x] = psdCompletion(x,K,cliqueDomain);

xVect = x(2:K.s(1),:);

WMat = reshape(x,K.s(1),K.s(1));

WMat = WMat(2:K.s(1),2:K.s(1)); 

return

