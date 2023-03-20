function [statusSW,costMatrix,maxDegree,minDegree,aveDegree] = ... 
    generateNetwork4(type,noOfVertices,nodeDegree,randSeed);

% type
%   = 1 ---> [0,1]*[0,1)
%   = 2 ---> 2-dimensional sphere
%   = 3 ---> torus

noOfVertices0 = noOfVertices;
noOfVertices = noOfVertices + 5*nodeDegree; 

% Generating random vertices in [0,1]x[0,1]
rand('state',randSeed);

if type == 1
% [0,1]x[0,1]
    xMatrix0 = rand(2,noOfVertices); 
    % Sorting the vertices according to the first coordinate
    [temp,perm] = sort(xMatrix0(1,:),'ascend'); 
    xMatrix0 = xMatrix0(:,perm);
elseif type == 2
% sphere
    xMatrix0 = (rand(3,noOfVertices) - 0.5);
    temp = xMatrix0 .* xMatrix0;
    xNorm = sqrt(sum(temp,1));
%     size(xNorm)
%     size(xMatrix0)
    xMatrix0 = xMatrix0 ./ (ones(3,1)*xNorm);
    [temp,perm] = sort(xMatrix0(1,:),'ascend'); 
    xMatrix0 = xMatrix0(:,perm);    
elseif type == 3
% torus
    xMatrix0 = (rand(2,noOfVertices0) - 0.5);
    temp = xMatrix0 .* xMatrix0;
    xNorm = sqrt(sum(temp,1));
    xMatrix0 = xMatrix0 ./ (ones(2,1)*xNorm);    
    xMatrix1 = (rand(2,noOfVertices0) - 0.5);
    temp = xMatrix1 .* xMatrix1;
    xNorm = sqrt(sum(temp,1));
    xMatrix1 = xMatrix1 ./ (ones(2,1)*xNorm);
    xMatrix1 = 0.2*xMatrix1; 
    xMatrix0 = [xMatrix0; xMatrix1];
    clear xMatrix1
    idxPlus = find(xMatrix0(2,:) > 0);
    idxMinus = find(xMatrix0(2,:) <= 0);
    xMatrixP = xMatrix0(:,idxPlus);
    [temp,perm] = sort(xMatrixP(1,:),'ascend'); 
    xMatrixP = xMatrixP(:,perm); 
    xMatrixM = xMatrix0(:,idxMinus);
    [temp,perm] = sort(xMatrixM(1,:),'descend'); 
    xMatrixM = xMatrixM(:,perm); 
%    xMatrix0 = [xMatrix0,xMatrix0(:,1:1+noOfVertices-noOfVertices0)];
    xMatrix0 = [xMatrixP,xMatrixM];
    noOfVertices = size(xMatrix0,2); 
else
    error('## type = 1, 2 or 3');
end

if type <= 2
    costMatrix = sparse(noOfVertices,noOfVertices);
    for q = 1:noOfVertices0
        oneRowIdx = repmat(q,1,noOfVertices);
        oneRowMat = xMatrix0-xMatrix0(:,oneRowIdx);
        oneRowMat = oneRowMat.*oneRowMat;
        %    distanceMatrix(q,q+1:noOfVertices) = sqrt(sum(oneRowMat,1));
        oneRow = sqrt(sum(oneRowMat,1));
        [temp,perm] = sort(oneRow(q+1:noOfVertices),'ascend');
        idx0 = perm(1:nodeDegree);
        idx1 = idx0+q;
        costMatrix(q,idx1) = oneRow(idx0);
        if q > 1
            [temp,idx] = min(oneRow(1:q-1));
            costMatrix(idx,q) = oneRow(idx);
        end
    end
    noOfVertices = noOfVertices0;
    xMatrix0 = xMatrix0(:,1:noOfVertices);
    costMatrix = costMatrix(1:noOfVertices,1:noOfVertices);
elseif type == 3
    nodeDegreePlus = nodeDegree;
    nodeDegreeMinus = 1; 
    costMatrix = sparse(noOfVertices,noOfVertices);
    nVhalf = floor(noOfVertices/2); 
    for q = 1:noOfVertices
        oneRowIdx = repmat(q,1,noOfVertices);
        oneRowMat = xMatrix0-xMatrix0(:,oneRowIdx);
        oneRowMat = oneRowMat.*oneRowMat;
        %    distanceMatrix(q,q+1:noOfVertices) = sqrt(sum(oneRowMat,1));
        oneRow = sqrt(sum(oneRowMat,1));
        [temp,perm] = sort(oneRow,'ascend');
        if q+nVhalf <= noOfVertices
            idxPlus = find((q < perm) & (perm <= q+nVhalf)); 
%            idxPlus = find((q < perm)) 
%            XXXX
            permPlus = perm(idxPlus); 
            idxMinus = find((perm < q) | (q+nVhalf < perm)); 
            permMinus = perm(idxMinus); 
        else
            idxMinus = find((perm < q) & (q-nVhalf < perm)); 
            permMinus = perm(idxMinus); 
            idxPlus = find((q < perm) | (perm <= q-nVhalf)); 
            permPlus = perm(idxPlus); 
        end
%         perm
%         nodeDegreePlus
        idxPlus = permPlus(1:nodeDegreePlus);   
        idxMinus = permMinus(1:nodeDegreeMinus);   
        costMatrix(q,idxPlus) = oneRow(idxPlus);
        costMatrix(q,idxMinus) = oneRow(idxMinus);
    end
    costMatrix = costMatrix + costMatrix';
    costMatrix = triu(costMatrix); 
end

% full(xMatrix0)
% full(costMatrix)

debugSW = 1;
if (debugSW == 1) && (type == 1)
    figNo = 11;
    figure(figNo);
    axis([0.0 1.0 0.0 1.0]);
    plot(xMatrix0(1,1:noOfVertices),xMatrix0(2,1:noOfVertices),'r.');
    hold on;
    fromNodes = [];
    toNodes = [];
    for q = 1:noOfVertices
        idx = find(costMatrix(q,:));
        r = length(idx);
        fromNodes = [fromNodes; repmat(xMatrix0(:,q)',r,1)];
        toNodes = [toNodes;xMatrix0(:,idx)'];
    end
    for i=1:size(fromNodes,1)
        xx = [fromNodes(i,1),toNodes(i,1)];
        yy = [fromNodes(i,2),toNodes(i,2)];
        plot(xx,yy,'b-');
    end
    hold off;
end

% full(costMatrix)

costMatrix = spones(costMatrix); 

% full(costMatrix)

nodeDegreeVector = full(sum(costMatrix+costMatrix',1));
maxDegree = max(nodeDegreeVector); 
minDegree = min(nodeDegreeVector);
aveDegree = sum(nodeDegreeVector)/length(nodeDegreeVector);

debugSW = 0;
if debugSW == 1
    sPattern = costMatrix + costMatrix' + (1+noOfVertices)*speye(noOfVertices,noOfVertices); 
    perm = symamd(sPattern);
    figure(type);
    spy(sPattern(perm,perm));    
%    spy(chol(sPattern(perm,perm)));    
end

statusSW = 1; 

return

    
