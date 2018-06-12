function demo_persistent_clique_conductance()

    setup;
    
    warning off
    data = load('data/HCP1200_Yeo100');
    A = data.ggmobj.stability_graphs(:,:,5);
    Ci = data.communityno;
    figure;
    imagesc(A); 
    axis image off;
    if(exist('brewermap'))
         colormap(flipud(brewermap([],'RdYlBu')));
    end
    colorbar;
    title('Schaefer Yeo 100 Network on HCP1200 Group Data'); 
    warning on;
    

    [metrics clique_adj cliques] = ...
         tda.persistent_conductance(A,Ci,false,'persistent_conductance');
    
    % Plot DMN-Other Communities Conductance
    create_matrix_movie(squeeze(metrics.conductances(:,:,:,7)), ...
                        '~/Downloads/HCP1200_DMN_Conductances'); 

    create_matrix_movie(squeeze(metrics.conductances(:,:,1:6,6)), ...
                             '~/Downloads/HCP1200_FPN_Conductances'); 
    