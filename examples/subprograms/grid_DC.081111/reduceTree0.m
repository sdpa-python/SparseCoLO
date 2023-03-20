%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clique,adjacencyMatrixT,edgeCostVectT,incidenceMatrixT] = ... 
    reduceTree0(clique,adjacencyMatrixT,edgeCostVectT,incidenceMatrixT,bdUnionOfCliques); 

noOfNodes = size(incidenceMatrixT,1); 
noOfEdges = size(incidenceMatrixT,2); 
nodeDegree = sum(spones(incidenceMatrixT),2); 

%%%%%%%%%%
debugSW =0;
if debugSW == 1
	clique
    fprintf('\nInformation on cliques: \n');        
    for p=1:clique.NoC
        fprintf(' %2d: ',p);
        for i=1:length(clique.Set{p})
            fprintf(' %2d',clique.Set{p}(i));
        end
        fprintf('\n');        
    end
    fprintf('Information on the clique tree computed\n');    
    fprintf('%d x %d adjacencyMatrixT\n',clique.NoC,clique.NoC)
    fprintf('     ');
    for q=1:clique.NoC
        fprintf('%3d',q);
    end
    fprintf('\n');
    for p=1:clique.NoC
        fprintf('%3d: ',p); 
        for q=1:clique.NoC
            fprintf('%3d',full(adjacencyMatrixT(p,q)));
        end
        fprintf('\n');
    end
    fprintf('           1 x %d edgeCostVectT\n',clique.NoC-1)
    fprintf('         ');
    for q=1:clique.NoC-1
        fprintf('%3d',full(edgeCostVectT(q)));
    end
    fprintf('\n');
    fprintf('     nDeg, %d x %d incidenceMatrixT\n',clique.NoC,clique.NoC-1)
    fprintf('         ');
    for q=1:clique.NoC-1
        fprintf('%3d',q);
    end
    fprintf('\n');
    for p=1:clique.NoC
        fprintf('%3d: ',p);
        fprintf('%3d ',nodeDegree(p)); 
        for q=1:clique.NoC-1
            fprintf('%3d',full(incidenceMatrixT(p,q)));
        end
        fprintf('\n');
    end
end
%%%%%%%%%

continueSW = 1;
noOfCliques = clique.NoC; 
pointer = 0; 
edgeIdxChecked = sparse(1,clique.NoC-1); 
combinePair = sparse(noOfCliques,1);
largeCliqueIdx = find(clique.NoElem >= bdUnionOfCliques); 
combinePair(largeCliqueIdx) = 1;
while continueSW == 1
    edgeIdx = find(edgeIdxChecked == 0);
    continueSW = 0; 
    if ~isempty(edgeIdx)
        for e=edgeIdx
            nodeIdx = find(incidenceMatrixT(:,e)');
            %            combinePair(nodeIdx,1)'
            %             e
            %             nodeIdx

            if (nnz(combinePair(nodeIdx,1)') == 0)
                tempCliqueSet = union(clique.Set{nodeIdx(1)},clique.Set{nodeIdx(2)});
                lenTempCliqueSet = length(tempCliqueSet);
                %
                %                 lenTempCliqueSet
                %                 bdUnionOfCliques

                if lenTempCliqueSet <= bdUnionOfCliques
                    continueSW = 1;
                    pointer = pointer + 1;
                    if lenTempCliqueSet >= bdUnionOfCliques
                        combinePair(nodeIdx(1),1) = 1;
                    end
                    combinePair(nodeIdx(2),1) = -1;
                    %                     format long
                    %                     full(combinePair)
                    %                     format short
                    edgeCostVectT(e) = -pointer;
                    clique.Set{nodeIdx(1)} = tempCliqueSet;
                    clique.NoElem(nodeIdx(1)) = lenTempCliqueSet;
                    clique.Set{nodeIdx(2)} = [];
                    clique.NoElem(nodeIdx(2)) = 0;
                    incidenceMatrixT(nodeIdx(1),:) = incidenceMatrixT(nodeIdx(1),:)+incidenceMatrixT(nodeIdx(2),:);
                    edgeIdxChecked(e) = -1;
                else
                    edgeIdxChecked(e) = 1;
                end
            end
            debugSW = 0;
            if debugSW == 1
                fprintf('\n');
                fprintf('edgeIdxChecked\n         ') 
                for q=1:clique.NoC-1
                	fprintf('%3d',full(edgeIdxChecked(q)));
                end
                fprintf('\n');
                fprintf('    cPair, %d x %d incidenceMatrixT\n',clique.NoC,clique.NoC-1)
                fprintf('         ');
                for q=1:clique.NoC-1
                    fprintf('%3d',q);
                end
                fprintf('\n');
                for p=1:clique.NoC
                    fprintf('%3d: ',p);
                    fprintf('%3d ',combinePair(p,1));
                    for q=1:clique.NoC-1
                        fprintf('%3d',full(incidenceMatrixT(p,q)));
                    end
                    fprintf('\n');
                end
            end
        end
%         newNodeSetToCheck = intersect(find(combinePair >= 0),find(combinePair <= 1.0e6)); 
%         combinePair(newNodeSetToCheck,1) = 0;
    end    
end

%XXXXX

% full(edgeIdxChecked)
% % 
% full(combinePair')


newEdgeIdx = find(edgeIdxChecked >= 0); 
newNodeIdx = find(combinePair' >= 0);

% edgeCostVectT
% full(incidenceMatrixT)

edgeCostVectT = edgeCostVectT(newEdgeIdx); 
incidenceMatrixT = incidenceMatrixT(newNodeIdx,newEdgeIdx); 

% edgeCostVectT 
% 
% full(incidenceMatrixT)

pointer = 0;
clique.Elem = [];
clique.NoElem = [];
for p=1:clique.NoC
    if ~isempty(clique.Set{p}) 
        pointer = pointer + 1;
        clique.Set{pointer} = clique.Set{p};
        clique.Elem = [clique.Elem,clique.Set{p}];
        clique.NoElem = [clique.NoElem, length(clique.Set{p})]; 
    end
end
for p=pointer+1:clique.NoC
    clique.Set{p} = [];
    clear clique.Set{p}
end
clique.NoC = size(incidenceMatrixT,1); 


clique.maxC = max(clique.NoElem);
clique.minC = min(clique.NoElem);

adjacencyMatrixT=sparse(clique.NoC,clique.NoC); 
for p=1:clique.NoC-1
    idx = find(incidenceMatrixT(:,p)'); 
    adjacencyMatrixT(idx(1),idx(2)) = edgeCostVectT(p);
end
nodeDegree = sum(spones(incidenceMatrixT),2); 

%%%%%%%%%%
debugSW =0;
if debugSW == 1
    clique
%     fprintf('\nInformation on cliques: \n');        
%     for p=1:clique.NoC
%         fprintf(' %2d: ',p);
%         for i=1:length(clique.Set{p})
%             fprintf(' %2d',clique.Set{p}(i));
%         end
%         fprintf('\n');        
%     end
    fprintf('Information on the clique tree computed\n');    
    fprintf('%d x %d adjacencyMatrixT\n',clique.NoC,clique.NoC)
    fprintf('     ');
    for q=1:clique.NoC
        fprintf('%3d',q);
    end
    fprintf('\n');
    for p=1:clique.NoC
        fprintf('%3d: ',p); 
        for q=1:clique.NoC
            fprintf('%3d',full(adjacencyMatrixT(p,q)));
        end
        fprintf('\n');
    end
    fprintf('           1 x %d edgeCostVectT\n',clique.NoC-1)
    fprintf('         ');
    for q=1:clique.NoC-1
        fprintf('%3d',full(edgeCostVectT(q)));
    end
    fprintf('\n');
    fprintf('     nDeg, %d x %d incidenceMatrixT\n',clique.NoC,clique.NoC-1)
    fprintf('         ');
    for q=1:clique.NoC-1
        fprintf('%3d',q);
    end
    fprintf('\n');
    for p=1:clique.NoC
        fprintf('%3d: ',p);
        fprintf('%3d ',nodeDegree(p)); 
        for q=1:clique.NoC-1
            fprintf('%3d',full(incidenceMatrixT(p,q)));
        end
        fprintf('\n');
    end
    XXXXX
end
%%%%%%%%%

% XXXXX

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
