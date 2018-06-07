function CC = clique_adjacency(cliques,n_nodes,k)
% Returns undirected clique motif co-occurances
% 
% INPUT
%  - cliques: cell array of cliques or n_cliques x n_nodes matrix
%  - n_nodes: number of nodes in network
%  - k : size of the clique motif (starts at 3 to k)
% OUTPUT
%   - CC: n_node x n_node x k-1. If CC(i,j,2) = 1 then i,j there is one 3-clique containing nodes i and j.
% 
% Example: 
% % Given A, a binary adjacency matrix
% [cliques] = all_cliques(A); 
% clique_adj = clique_adjacency(cliques,length(A),3);

    
    cellmode = iscell(cliques);
    if(cellmode)
        CC = zeros(n_nodes,n_nodes,k-1);
        n_cliques = length(cliques);
    else
        CC = zeros(n_nodes,n_nodes,k-1);
        [n_cliques n_nodes] = size(cliques);
        clique_size = sum(cliques,2);        
    end

    for clique_no=1:n_cliques
        for kk=2:k 
            if(cellmode && kk>=3)       
                if(length(cliques{clique_no})==kk)
                    tmp_CC = zeros(n_nodes,n_nodes);
                    tmp_CC(cliques{clique_no},cliques{clique_no}) = 1;
                    CC(:,:,kk-1) = CC(:,:,kk-1) + tmp_CC;
                end
            else
                if(clique_size(clique_no)==kk)
                    tmp_CC = zeros(n_nodes,n_nodes);
                    clique_nodes = find(cliques(clique_no,:));
                    tmp_CC(clique_nodes,clique_nodes) = 1;
                    CC(:,:,kk-1) = CC(:,:,kk-1) + tmp_CC;
                end
            end
        end
    end
    
end