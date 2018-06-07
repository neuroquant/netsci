function [metrics clique_adj clique_array] = persistent_conductance(A,Ci,varargin)
    % Conductance between source and target communities s,t using clique adjacency information
    % 
    % USAGE
    %   [metrics clique_adj cliques] = persistent_conductance(Sig,Ci)
    % INPUT
    %   A  - A weighted or binary adjacency matrix of p x p nodes
    %   Ci - A p x 1 vector of community partition labels for each node
    % 
    % OUTPUT
    %
    % metrics
    %   metrics.conductances 
    %   metrics.degree
    % 
    % clique_adj - each cell in array contains node x node x clique_size adjacency
    %
    % cliques   - each cell in array contains cliques x node members
    import +community.*
    
    addpath('~/MATLAB/netsci/external/k_clique');
    addpath('~/MATLAB/packages/matlab-bgl');
    
    %% For stability graphs
    % thresholds = fliplr(round(linspace(.5,.8,20),2));
    % For (-1,1) similarity or correlations
    
    edges = A(find(triu(A,1)));
    edge_range = prctile(abs(edges),[75 99]);
    min_edge = edge_range(1); 
    max_edge = edge_range(2);     
    thresholds = fliplr(round(linspace(min_edge,max_edge,50),2));
    
    loaddata = false;
    use_clique_top = true;
    k_motifs = 15; % starts from 2 to k_motifs+1
    n_nodes = length(A);
    communities = unique(Ci);
    n_communities = length(communities); 
    conductances = ...
     zeros(length(thresholds),k_motifs,n_communities,n_communities);
    degree = zeros(length(thresholds),k_motifs,n_nodes);
    volume = zeros(length(thresholds),k_motifs,n_nodes,n_communities); 
    mstree = zeros(length(thresholds),n_nodes);
    components = zeros(length(thresholds),n_nodes,10);
    betti = zeros(length(thresholds),2);
    births = zeros(length(thresholds),k_motifs,n_nodes,n_nodes);
    
    A = (A + A')/2;
    for th_no=1:length(thresholds)
        th_no
        if(loaddata)
            load('tmp/persistant_conductance', ...
              'clique_array','clique_adj');
        else 
            if(use_clique_top)
                
                clique_array{th_no} = ...
                            community.cuts.clique_top(A,thresholds(th_no));
                clique_adj{th_no} = ...
                            clique_adjacency(clique_array{th_no},...
                                                n_nodes, ...
                                                k_motifs+1 ...
                                                );
                
            else
                clique_array{th_no} = ...
                     cuts.get_cliques(A,thresholds(th_no));
                clique_adj{th_no} = clique_adjacency(clique_array{th_no}, ...
                                            n_nodes,k_motifs+1);
            end
                                            
            clique_adj{th_no}(:,:,1) = (abs(A)>=thresholds(th_no));
            
        end
        % conductances(th_no,1) = ...
%             cuts.conductance(A.*(abs(A)>=thresholds(th_no)),Ci==6,Ci==7);
        for kk=1:k_motifs
            degree(th_no,kk,:) = sum(mean(clique_adj{th_no}(:,:,1:kk)~=0),3);
            for Cii=1:n_communities
                volume(th_no,kk,:,Cii) = community.cuts.volume(...
                            mean(clique_adj{th_no}(:,:,kk),3), ...
                            Ci==Cii);
                for Cjj=Cii:n_communities
                    conductances(th_no,kk,Cii,Cjj) = ...
                             community.cuts.conductance( ...
                                mean(clique_adj{th_no}(:,:,kk),3), ...
                                Ci==Cii, ...
                                Ci==Cjj);
                end
            end
        end
        tmp_A = (clique_adj{th_no}(:,:,1).*(A));
        %tmp_mst = mst(tmp_A);
        [H T] = dendrogram(linkage(tmp_A));
        mstree(th_no,:) = T;
        for ii=1:max(T)
            component_sz(ii) = sum(T==ii);
        end
        for ii=1:size(betti,3)
            components(th_no,:,ii) = (T==ii);
        end
        betti(th_no,1) = max(T);
        betti(th_no,2) = sum(sum(clique_adj{th_no}(:,:,1))) ...
                                 - n_nodes + betti(th_no,1); 
    end
    

    metrics = {};
    metrics.conductances = conductances;
    metrics.degree = degree;
    metrics.volume = volume;
    metrics.mstree = mstree;
    metrics.betti = betti; 
    
    if(~exist('tmp'))
        mkdir('tmp');
    end
    
    save('tmp/persistant_conductance', ...
      'metrics','clique_array','clique_adj','k_motifs','A','Ci');    
    
end

function opts = create_options(varargin)
    
    
    
end

