function [] = plot_behavior_under_partition(x, labels, CVI, CriterionValue_inc, nClusters_inc, valind_inc, ART, x_axis_ticks, CVI_name, class_order, classes_merged, y)

    %% Setup 
    [~, dim] = size(x);  
    experiment = 'Under-partitioned data';
    
    % Extract relevant data
    [x_samples, y_cvi, y_clusters, x_after, offset, start, stop] = extract_under_partition_data(CriterionValue_inc, nClusters_inc, x_axis_ticks, class_order, classes_merged, y);  
    
    % Under-partition time interval
    UP_cluster_start = start;
    UP_cluster_stop = stop;

    % Event position
    under_position = x_after(1) + offset;

    % Ticks  
    inds1 = x_axis_ticks < under_position;
    ind2 = find(inds1, 1, 'last');
    XTick = [x_axis_ticks(1:ind2) under_position x_axis_ticks(ind2+1:end)]; 
    XTickLabels = num2cell(XTick);
    XTickLabels{ind2} = ''; 
    XTickLabels{ind2+1} = 'UP'; 
        
    % Figure parameters
    FONTSIZE = 16;
    FONTWEIGHT = 'bold';
    LINEWIDTH = 2;
    
    % Open figure
    figure('visible', 'on');    
    set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) % Fullscreen

    %% iCVI 
    subplot(2,2,1);     
    box on
    hold on
    title(CVI_name, 'Interpreter', 'none')
    maximum = max(y_cvi);
    minimum = min(y_cvi);
    p1 = plot(x_samples, y_cvi, 'b', 'display', ['i' CVI], 'Linewidth', LINEWIDTH);
    if minimum ~= maximum
        ylim([minimum maximum])
    end            
    xlim([min(x_samples)-1 max(x_samples)+1]) 
    xlabel('Time', 'FontWeight', FONTWEIGHT)
    ylabel('iCVI', 'FontWeight', FONTWEIGHT)     

    % Axis properties
    box on
    set(gcf, 'Color', 'w'); 
    set(gca, 'Color', 'none'); 
    ax = gca;
    ax.GridLineStyle = '--';
    ax.GridColor = [0 0 0];
    ax.GridAlpha=.5;
    ax.FontSize = FONTSIZE;
    ax.FontWeight = FONTWEIGHT;
    ax.XGrid = 'on';
    xtickangle(ax, 90);
    set(gca, 'XTick', XTick)  
    set(gca, 'XTickLabels', XTickLabels)    
    
    % Vertical line on event position
    xval = [under_position, under_position];    
    yval = [minimum, maximum];    
    plot(xval, yval, 'Color', [0 .8 0], 'Linewidth', 2, 'LineStyle', '-')

    % Patch
    min_y = minimum;
    max_y = maximum;
    x_axis = [UP_cluster_start UP_cluster_stop UP_cluster_stop UP_cluster_start];
    y_axis = [min_y min_y max_y max_y];
    p2 = patch(x_axis, y_axis, 'black');
    p2.EdgeColor = 'none';
    p2.FaceAlpha = 0.2;

    % Order of plots
    uistack(p1,'top');

    %% Number of Clusters   
    subplot(2,2,3);     
    box on
    hold on
    maximum = max(y_clusters);
    minimum = min(y_clusters);
    p1 = plot(x_samples, y_clusters, 'r', 'Linewidth', LINEWIDTH);
    ylim([minimum-1 maximum+1])
    xlim([min(x_samples)-1 max(x_samples)+1]) 
    ylabel('No. Clusters','FontWeight',FONTWEIGHT)
    xlabel('Time','FontWeight',FONTWEIGHT)
         
    % Axis properties
    box on
    set(gcf, 'Color','w'); 
    set(gca, 'Color', 'none');
    ax = gca;
    ax.GridLineStyle = '--';
    ax.GridColor = [0 0 0];
    ax.GridAlpha=.5;
    ax.FontSize = FONTSIZE;
    ax.FontWeight = FONTWEIGHT;
    ax.XGrid = 'on';
    xtickangle(ax, 90);
    set(gca, 'XTick', XTick)  
    set(gca, 'XTickLabels', XTickLabels) 
    
    % Vertical line on event position
    xval = [under_position, under_position];    
    yval = [minimum-1, maximum+1];   
    plot(xval,yval,'Color',[0 .8 0],'Linewidth', 2, 'LineStyle', '-')

    % Patch
    min_y = minimum-1;
    max_y = maximum+1;
    x_axis = [UP_cluster_start UP_cluster_stop UP_cluster_stop UP_cluster_start];
    y_axis = [min_y min_y max_y max_y];
    p2 = patch(x_axis, y_axis, 'black');
    p2.EdgeColor = 'none';
    p2.FaceAlpha = 0.2;

    % Order of plots
    uistack(p1,'top');
    
    %% Data
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
    set(gcf, 'Color', 'w'); 
    set(gca, 'Color', 'none'); % sets axes background
    ax = gca;
    ax.FontSize = FONTSIZE;
    ax.FontWeight = FONTWEIGHT;   

end