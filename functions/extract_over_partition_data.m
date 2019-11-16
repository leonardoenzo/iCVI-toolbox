function [x_samples, y_cvi, y_clusters, x_tick_labels, x_after, offset] = extract_over_partition_data(CriterionValue_inc, nClusters_inc, x_axis_ticks, class_selected, y)

    nSamples = length(nClusters_inc);

    % Remove 1st cluster
    ind = nClusters_inc<=1;
    offset = sum(ind);
    x_samples = 1:nSamples;
    x_samples(ind) = [];
    y_cvi = CriterionValue_inc;
    y_cvi(ind) = [];
    y_clusters = nClusters_inc;
    y_clusters(ind) = [];
    
    % Get xticks
    x_tick_labels = x_axis_ticks;
    
    % Search splitting time
    start = find(y==class_selected, 1, 'first');
    stop = find(y==class_selected, 1, 'last');        
    for tx=start:stop-1                      
        if nClusters_inc(tx+1) > nClusters_inc(tx)
            t = [tx tx+1];
        end
    end   
    
    % Set before/after times
    x_before(1) = start - offset;
    x_before(2) = t(1) - offset;
    x_after(1)  = t(2) - offset;
    x_after(2)  = stop - offset;   
    
end