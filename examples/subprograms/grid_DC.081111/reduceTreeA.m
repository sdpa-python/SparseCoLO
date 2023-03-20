%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [modifiedSW,clique,adjacencyMatrixT,edgeCostVectT,incidenceMatrixT] = ... 
    reduceTreeA(clique,adjacencyMatrixT,edgeCostVectT,incidenceMatrixT,sigma0); 

modifiedSW = 0; 
noOfNodes = size(incidenceMatrixT,1); 
noOfEdges = size(incidenceMatrixT,2); 
nodeDegree = sum(spones(incidenceMatrixT),2); 
activeNode = [1:noOfNodes];
activeEdge = [1:noOfEdges]; 

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

% full(nodeDegree)

idxNode3 = find(nodeDegree' >= 3); 
for i1=idxNode3
    idxEdge3 = find(incidenceMatrixT(i1,:));    
    for j=idxEdge3
%         j
%         size(incidenceMatrixT)
%         
        adjTwoNodes = find(incidenceMatrixT(:,j)');       
        i2 = adjTwoNodes(1); 
        if i1 == i2
            i2 = adjTwoNodes(2); 
        end
        %%%%%%%%%%
        debugSW = 0; 
        if (debugSW == 1) && (i1==10) && (j==7)
            idxEdge3
            i1
            j
            adjTwoNodes
            i2
            nodeDegree(i2)
%             clique
%             clique.Elem
%             for p=1:clique.NoC
%                 clique.Set{p}
%             end
%            XXXXX
        end
        %%%%%%%%%
        sigma = min([edgeCostVectT(j)/clique.NoElem(i1),edgeCostVectT(j)/clique.NoElem(i2)]); 
        if (nodeDegree(i2) == 1) && (sigma > sigma0) 
            modifiedSW = 1; 
            incidenceMatrixT(i1,:) = incidenceMatrixT(i1,:) + incidenceMatrixT(i2,:); 
            incidenceMatrixT(:,j) = sparse(noOfNodes,1); 
            incidenceMatrixT(i2,:) = sparse(1,noOfEdges); 
            activeEdge = setdiff(activeEdge,j);
            activeNode = setdiff(activeNode,i2);
            clique.Set{i1} = union(clique.Set{i1},clique.Set{i2});
            clique.Set{i1} = sort(clique.Set{i1});
            clique.NoElem(i1) = length(clique.Set{i1}); 
            clique.Set{i2} = []; 
            clique.NoElem(i2) = 0;
%             nonemptyClIdx = find(clique.NoElem); 
%             clique.maxC = max(clique.maxC,clique.NoElem(i1)); 
%             clique.NoElem = clique.NoElem(nonemptyClIdx); 
%             clique.minC = min(clique.NoElem); 
%             pointer = 0;
%             clique.Elem = [];
%             for p=1:clique.NoC
%                 if p ~= i2
%                     pointer = pointer + 1;
%                     clique.Set{pointer} = clique.Set{p};
%                     clique.Elem = [clique.Elem,clique.Set{p}];
%                 end
%             end
%             clique.Set{clique.NoC} = [];
%             clear clique.Set{clique.NoC}; 
%            clique.NoC = clique.NoC - 1;
            debugSW = 0;
            if (debugSW == 1) && (i1==10) && (j==7) 
%                 clique
%                 clique.Elem
                for p=1:clique.NoC
                    clique.Set{p}
                end
                activeEdge
                full(incidenceMatrixT)
                activeNode
                XXXX
            end
        end
    end
end
incidenceMatrixT = incidenceMatrixT(activeNode,activeEdge); 
edgeCostVectT = edgeCostVectT(1,activeEdge); 
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
clique.NoC = length(activeNode); 
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
