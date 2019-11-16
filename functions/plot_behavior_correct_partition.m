function [] = plot_behavior_correct_partition(x, labels, CVI, CriterionValue_inc, nClusters_inc, valind_inc, ART, x_axis_ticks, CVI_name)

    [nSamples, dim] = size(x);   
    experiment = 'Correctly partitioned data';

    % Figure parameters
    FONTSIZE = 16;
    FONTWEIGHT = 'bold';
    LINEWIDTH = 2;

    % Plot new figure
    figure('visible', 'on');
    
    % Remove 1st cluster
    ind = nClusters_inc<=1;
    var_x = 1:nSamples;
    var_x(ind) = [];
    var_y1 = CriterionValue_inc;
    var_y1(ind) = [];
    var_y2 = nClusters_inc;
    var_y2(ind) = [];
    
    % CVI
    subplot(2,2,1);
    title(CVI_name, 'Interpreter', 'none')
    hold on   
              
    maximum = max(var_y1);
    minimum = min(var_y1);     
    plot(var_x,var_y1,'b', 'display', ['i' CVI], 'Linewidth', LINEWIDTH);              
    if minimum ~= maximum
        ylim([minimum maximum])
    end                  
    xlabel('Time','FontWeight', FONTWEIGHT)
    ylabel('iCVI','FontWeight',FONTWEIGHT)   

    xlim([min(var_x) max(var_x)]) 

    % Axis properties
    box on
    set(gcf,'color','w'); 
    set(gca, 'XTick', x_axis_ticks)
    set(gca, 'Color', 'none'); 
    ax = gca;
    ax.GridLineStyle = '--';
    ax.GridColor = [0 0 0];
    ax.GridAlpha=.5;
    ax.FontSize = FONTSIZE;
    ax.FontWeight = FONTWEIGHT;
    ax.XGrid = 'on';
    xtickangle(ax, 90);
    
    % Number of clusters
    subplot(2,2,3);
    hold on   
             
    plot(var_x, var_y2, 'r', 'Linewidth', LINEWIDTH); 
    ylabel('No. Clusters', 'FontWeight', FONTWEIGHT)
    xlabel('Time','FontWeight', FONTWEIGHT)  

    ylim([1 max(var_y2)+1])
    xlim([min(var_x) max(var_x)]) 

    % Axis properties
    box on
    set(gcf,'color','w'); 
    set(gca, 'XTick', x_axis_ticks)
    set(gca, 'Color', 'none'); 
    ax = gca;
    ax.GridLineStyle = '--';
    ax.GridColor = [0 0 0];
    ax.GridAlpha=.5;
    ax.FontSize = FONTSIZE;
    ax.FontWeight = FONTWEIGHT;
    ax.XGrid = 'on';
    xtickangle(ax, 90);
    
    % Data
    subplot(1,2,2)
    title(experiment)
    hold on   
    if dim==2 
        clrs = rand(nClusters_inc(end), 3);
        if strcmp(CVI, 'CONN')    
            % Prototypes & data for CONN_Index
            draw_FA_categories(ART.ART_A, x, labels, 2, ART.map, false, clrs)
            hold on
            CONNvis(ART.ART_A, ART.map, valind_inc.CONN, clrs);
            set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]);  
        else            
            for ix=1:nClusters_inc(end)
                plot(x(labels==ix,1),x(labels==ix,2), 'Marker', 'o', 'MarkerSize', 10, 'Color', clrs(ix,:), 'LineStyle', 'None');
            end
        end
    end
    % Axis properties
    box on
    set(gcf, 'Color','w'); 
    set(gca, 'Color', 'none');
    ax = gca;
    ax.FontSize = FONTSIZE;
    ax.FontWeight = FONTWEIGHT;
    
    % Fullscreen
    set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1])   
    
end