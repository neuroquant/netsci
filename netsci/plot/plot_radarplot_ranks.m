function figobj = plot_radarplot_ranks(Xranks,varargin)
    % PLOT_RADARPLOT_RANKS
    % 
    % - Xranks is a observation x feature matrix
    % 
    % 
    
    if(~exist('radar_plot'))
        error('Please add radarplot package')
    end
    
    [n p] = size(Xranks);
    maxrank = max(max(Xranks));
    Xranks = Xranks./(maxrank+1);
    
    if(nargin<2)
        %dim_legend = strcat({'Metric '},num2str([1:p]'));
        dim_legend = {'Bet-G','Bet-C','Bet-R', ...
                        'Close-G','Close-C','Close-R', ...
                        'Eff-G','Eff-C','Eff-R'};
        savefilename = 'tmp/radar_plots.pdf'                
    elseif(nargin<3)
        dim_legend = varargin{1};
        savefilename = 'tmp/radar_plots.pdf'
    else
        dim_legend = varargin{1};
        savefilename = [varargin{2} '.pdf'];
    end
    
    data_legend = repmat({'o'},[1 n])';
    params = {};
    %params.axes_handle = gca;
    %params.axes_visible = true;
    params.axis_lims = repmat([0 1],[p 1]);
    %params.axes_color = [0 1];
    %params.axis_reverse = [];
    %params.outer_isocurve_color = [];
    %params.num_isocurves = 5;
    params.circular_isocurves = true;
    params.datalines_transparency_value = .8;
    params.isocurve_line_style = '-'; %'-' | '--' | ':' | '-.' | 'none'
    params.show_dim_scale_label = false;
    params.rotate_dim_legend = true;
    %params.dim_legend_font = Helvetica;
    %params.dim_legend_font_bold = true;
    params.dim_legend_font_size = 14;
    params.data_symbols = data_legend;

    for ii=1:n
        h_tmp = figure('visible','off');
        % panelobj = panel(h_tmp);
        params.axes_handle = gca;
        figobj(ii) = radar_plot(Xranks(ii,:)',dim_legend,data_legend,params);
        if(ii==1)
            export_fig(savefilename,'-q100');
        else
            export_fig(savefilename,'-q100','-append');
        end
        close;
    end
    %panelobj.pack()
    
    %panelobj.add(figobj)
    
    
end