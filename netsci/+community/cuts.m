classdef cuts
% 
% Conductance is a quality function used in graph and spectral clustering
%
% REFERENCES
% 
%  J. Shi and J. Malik. Normalized cuts and image
% segmentation. IEEE Transcations of Pattern Analysis and
% Machine Intelligence, 22(8):888–905, 2000.
% 
% R. Kannan, S. Vempala, and A. Vetta. On clusterings:
% Good, bad and spectral. Jour. of the ACM, 51(3), 2004.
%
% Gleich, David. "Hierarchical directed spectral graph partitioning." Information Networks (2006).

	methods(Static)

		function output = conductance(A,S,T,varargin)
        %
        % INPUTS
        % - A is a n x n adjacency matrix (binary or weighted | undirected) 
		% - S is a n x 1 vector of node community membership or graph partitions. Community memberships indexed by 1,...,C or just a vector of 1s and 0s for a single group or cluster
        % - T is a n x 1 vector of 1s and 0s indicating cluster different from S. 
        
    	    switch nargin 
            case 4
                opts = varargin{1};
            case 5
                opts = varargin{1};
                warning('4th argument ignored');
            otherwise
                opts = community.cuts.create_options();
    		end
        
            nodelist = 1:length(A);
            n_partitions = max(S);
            % By default, let T be complement to 1st partition
            if(isempty(T) | n_partitions > 1)
                disp('Mode: All Pairwise Conductances')
                conductance = zeros(n_partitions,n_partitions);
                for cluster_i=1:n_partitions
                    T = 1*(S~=cluster_i);
                    vol_S = community.cuts.volume(A,1*(S==cluster_i));
                    vol_T = community.cuts.volume(A,T);
                    conductance(cluster_i,cluster_i) = ...
                                sum(sum(A(find(S),find(T))))/min(vol_S,vol_T); 
                    clear T vol_T;
                    for cluster_j=(cluster_i+1):n_partitions
                        T = 1*(S==cluster_j);
                        vol_T = community.cuts.volume(A,T);
                        conductance(cluster_i,cluster_j) = ...
                                sum(sum(A(find(S),find(T))))/min(vol_S,vol_T);
                        clear T vol_T;
                    end
                
                end
            else                
                vol_S = community.cuts.volume(A,S);
                vol_T = community.cuts.volume(A,T);                
                conductance = ...
                            sum(sum(A(find(S),find(T))))/min(vol_S,vol_T); 
            end
            output = conductance;
        end


		function output =  cut_size(A,S,T,varargin)
        %
        % INPUTS
        % - A is an adjacency matrix (binary or weighted | undirected) 
		% - S is a vector of node community membership or graph partitions


	
		end


		function output = volume(A,S,varargin)
        %
        % INPUTS
        % - A is an adjacency matrix (binary or weighted | undirected) 
		% - S is a 0/1 vector indicating cluster membership
        
            vol_nodes = sum(A,2);
            output = sum(vol_nodes(find(S)));
	
		end
		
        function options = create_options()
            
            options = {};
            options.isBinary = true; 
            options.pairwise = false; % if false then S vs. V \ S
            options.motif = 'edge'; % 2-clique, 3-clique, 4-cliques
             
        end
    
        function clique_mat = clique_top(A,thresh,varargin)
            
            find_max_cliques = true;
            min_clique = 2;
            max_clique = 20;
            
            [~, clique_mat] = ...
                 Cliquer.FindAll(1.0*(abs(A)>=thresh),...
                                 min_clique, ...
                                 max_clique, ...
                                 find_max_cliques, ...
                                 1000 ...
                                 );
            
            
        end
    
    
        function clique_array = get_cliques(A,thresh)
        
           tmpdir = fullfile('~/MATLAB/continuity_netsci','tmp')
           unix(['mkdir -p ' tmpdir]);
           
           G =  graph((abs(A)>=thresh));
           writetable(G.Edges, fullfile(tmpdir,'tmp_edgelist.txt'), ...
                      'WriteVariableNames',false,'Delimiter',' ');
        
           clique_no = 0; % all cliques computed regardless
           unix(['rm -r ' fullfile(tmpdir,'tmp_edgelist')])
           currdir = pwd;
           cd('~/MATLAB/netsci/external/CFinder-v2.0.6/')
           if(clique_no~=0)             
               clique_cmd = ...
                sprintf(['./CFinder ' ...
                ... %'-l ~/MATLAB/netsci/external/CFinder-v2.0.6/ ' ...
                '-i ~/MATLAB/continuity_netsci/tmp/tmp_edgelist.txt ' ...
                '-o ~/MATLAB/continuity_netsci/tmp/tmp_edgelist ' ...
                ' -k ' num2str(clique_no)]);
           else
                clique_cmd = ...
                 sprintf(['./CFinder ' ...
                 ... %'-l ~/MATLAB/netsci/external/CFinder-v2.0.6/ ' ...
                 '-i ~/MATLAB/continuity_netsci/tmp/tmp_edgelist.txt ' ...
                 '-o ~/MATLAB/continuity_netsci/tmp/tmp_edgelist ']);
            end
           unix(clique_cmd);
           cd(currdir)
           % Read CFinder Results
           cliquedir = '~/MATLAB/continuity_netsci/tmp/tmp_edgelist';
           cliquefile = fullfile(cliquedir,'cliques');
           % clique_array = readtable(cliquefile,...
           %      'ReadVariableNames',0, ...
           %      'HeaderLines',6,...
           %      'ReadRowNames',1,...
           %      'Delimiter',' ');
           try
               clique_mat = dlmread(cliquefile,' ',6,1);
           catch me
               disp(me)
               warning('File is empty')
           end
           clique_array = {};
           for ii=1:size(clique_mat,1)
               clique_array{ii} = setdiff(unique(clique_mat(ii,:)),0);
           end
        
        end
        
	end
end
	
	% def cut_size(G, S, T=None, weight=None):
	%     """Returns the size of the cut between two sets of nodes.
	%
	%     A *cut* is a partition of the nodes of a graph into two sets. The
	%     *cut size* is the sum of the weights of the edges "between" the two
	%     sets of nodes.
	%
	%     Parameters
	%     ----------
	%     G : NetworkX graph
	%
	%     S : sequence
	%         A sequence of nodes in `G`.
	%
	%     T : sequence
	%         A sequence of nodes in `G`. If not specified, this is taken to
	%         be the set complement of `S`.
	%
	%     weight : object
	%         Edge attribute key to use as weight. If not specified, edges
	%         have weight one.
	%
	%     Returns
	%     -------
	%     number
	%         Total weight of all edges from nodes in set `S` to nodes in
	%         set `T` (and, in the case of directed graphs, all edges from
	%         nodes in `T` to nodes in `S`).
	%
	%     Examples
	%     --------
	%     In the graph with two cliques joined by a single edges, the natural
	%     bipartition of the graph into two blocks, one for each clique,
	%     yields a cut of weight one::
	%
	%         >>> G = nx.barbell_graph(3, 0)
	%         >>> S = {0, 1, 2}
	%         >>> T = {3, 4, 5}
	%         >>> nx.cut_size(G, S, T)
	%         1
	%
	%     Each parallel edge in a multigraph is counted when determining the
	%     cut size::
	%
	%         >>> G = nx.MultiGraph(['ab', 'ab'])
	%         >>> S = {'a'}
	%         >>> T = {'b'}
	%         >>> nx.cut_size(G, S, T)
	%         2
	%
	%     Notes
	%     -----
	%     In a multigraph, the cut size is the total weight of edges including
	%     multiplicity.
	%
	%     """
	%     edges = nx.edge_boundary(G, S, T, data=weight, default=1)
	%     if G.is_directed():
	%         edges = chain(edges, nx.edge_boundary(G, T, S, data=weight, default=1))
	%     return sum(weight for u, v, weight in edges)
	%
	%
	%
	% [docs]def volume(G, S, weight=None):
	%     """Returns the volume of a set of nodes.
	%
	%     The *volume* of a set *S* is the sum of the (out-)degrees of nodes
	%     in *S* (taking into account parallel edges in multigraphs). [1]
	%
	%     Parameters
	%     ----------
	%     G : NetworkX graph
	%
	%     S : sequence
	%         A sequence of nodes in `G`.
	%
	%     weight : object
	%         Edge attribute key to use as weight. If not specified, edges
	%         have weight one.
	%
	%     Returns
	%     -------
	%     number
	%         The volume of the set of nodes represented by `S` in the graph
	%         `G`.
	%
	%     See also
	%     --------
	%     conductance
	%     cut_size
	%     edge_expansion
	%     edge_boundary
	%     normalized_cut_size
	%
	%     References
	%     ----------
	%     .. [1] David Gleich.
	%            *Hierarchical Directed Spectral Graph Partitioning*.
	%            <https://www.cs.purdue.edu/homes/dgleich/publications/Gleich%202005%20-%20hierarchical%20directed%20spectral.pdf>
	%
	%     """
	%     degree = G.out_degree if G.is_directed() else G.degree
	%     return sum(d for v, d in degree(S, weight=weight))
	%
	%
	%
	% [docs]def normalized_cut_size(G, S, T=None, weight=None):
	%     """Returns the normalized size of the cut between two sets of nodes.
	%
	%     The *normalized cut size* is the cut size times the sum of the
	%     reciprocal sizes of the volumes of the two sets. [1]
	%
	%     Parameters
	%     ----------
	%     G : NetworkX graph
	%
	%     S : sequence
	%         A sequence of nodes in `G`.
	%
	%     T : sequence
	%         A sequence of nodes in `G`.
	%
	%     weight : object
	%         Edge attribute key to use as weight. If not specified, edges
	%         have weight one.
	%
	%     Returns
	%     -------
	%     number
	%         The normalized cut size between the two sets `S` and `T`.
	%
	%     Notes
	%     -----
	%     In a multigraph, the cut size is the total weight of edges including
	%     multiplicity.
	%
	%     See also
	%     --------
	%     conductance
	%     cut_size
	%     edge_expansion
	%     volume
	%
	%     References
	%     ----------
	%     .. [1] David Gleich.
	%            *Hierarchical Directed Spectral Graph Partitioning*.
	%            <https://www.cs.purdue.edu/homes/dgleich/publications/Gleich%202005%20-%20hierarchical%20directed%20spectral.pdf>
	%
	%     """
	%     if T is None:
	%         T = set(G) - set(S)
	%     num_cut_edges = cut_size(G, S, T=T, weight=weight)
	%     volume_S = volume(G, S, weight=weight)
	%     volume_T = volume(G, T, weight=weight)
	%     return num_cut_edges * ((1 / volume_S) + (1 / volume_T))
	%
	%
	%
	% [docs]def conductance(G, S, T=None, weight=None):
	%     """Returns the conductance of two sets of nodes.
	%
	%     The *conductance* is the quotient of the cut size and the smaller of
	%     the volumes of the two sets. [1]
	%
	%     Parameters
	%     ----------
	%     G : NetworkX graph
	%
	%     S : sequence
	%         A sequence of nodes in `G`.
	%
	%     T : sequence
	%         A sequence of nodes in `G`.
	%
	%     weight : object
	%         Edge attribute key to use as weight. If not specified, edges
	%         have weight one.
	%
	%     Returns
	%     -------
	%     number
	%         The conductance between the two sets `S` and `T`.
	%
	%     See also
	%     --------
	%     cut_size
	%     edge_expansion
	%     normalized_cut_size
	%     volume
	%
	%     References
	%     ----------
	%     .. [1] David Gleich.
	%            *Hierarchical Directed Spectral Graph Partitioning*.
	%            <https://www.cs.purdue.edu/homes/dgleich/publications/Gleich%202005%20-%20hierarchical%20directed%20spectral.pdf>
	%
	%     """
	%     if T is None:
	%         T = set(G) - set(S)
	%     num_cut_edges = cut_size(G, S, T, weight=weight)
	%     volume_S = volume(G, S, weight=weight)
	%     volume_T = volume(G, T, weight=weight)
	%     return num_cut_edges / min(volume_S, volume_T)

	