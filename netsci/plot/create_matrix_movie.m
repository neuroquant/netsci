function create_matrix_movie(A,movie_filename)
    % 
    % A - p x p x nframes
    % movie_filename - where to save .avi
    
    if(isempty(movie_filename))
        movie_filename = fullfile('tmp','matrix_movie');
    end
    
    v = VideoWriter([movie_filename '.avi']);
    open(v);
    
    if(size(A,1)==size(A,2))
        imagesc(A(:,:,1))
    else
        imagesc(A(:,:,1)');
    end
    axis tight;
    colorbar;
    set(gcf,'Position',[60 75 450 350]);
    set(gca,'nextplot','replacechildren'); 
    
    %Create a set of frames and write each frame to the file.
    for k = 1:size(A,3)
       if(size(A,1)==size(A,2))
           imagesc(A(:,:,k));
           axis image off;
       else
           imagesc(A(:,:,k)');
           xlabel('Thresholds (Large -> Small)');
           ylabel('Clique Size');  
           title(['Frame No. ' num2str(k)]);
           axis tight;
           set(gca,'FontSize',18, ...
               'YTick',[2:2:15]','YTickLabel',num2str([3:2:16]'));
       end
       colormap(flipud(brewermap([],'RdYlBu')));
       colorbar;
       pause(.5);
       
       % ax = gca;
       % ax.Units = 'pixels';
       % pos = ax.Position;
       % marg = 30;
       % rect = [-marg+pos(1), -marg+pos(2), pos(3)+2*marg, pos(4)+2*marg];
       % ax.Units = 'normalized';
       frame = getframe(gcf);
       writeVideo(v,frame);
    end
    close(v);
end