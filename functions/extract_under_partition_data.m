function [x_samples, y_cvi, y_clusters, x_after, offset, start, stop] = extract_under_partition_data(CriterionValue_inc, nClusters_inc, x_axis_ticks, class_order, classes_merged, y)    

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
   
    % Search under-partition time
    ind1 = find(class_order == classes_merged(1));  % index of the merged class 1 in the suffled stream
    ind2 = find(class_order == classes_merged(2));  % index of the merged class 2 in the suffled stream       
    merge_time = class_order(max(ind1, ind2));      % get which merged class will be presented second (where the problem will be observable)
    start = find(y==merge_time, 1, 'first');
    stop = find(y==merge_time, 1, 'last');
    
    % Set before/after times
    ind3 = find(x_axis_ticks == start-1) - 1; % ind3: last sample of previous+previous cluster   
    x_before(1) = x_axis_ticks(ind3) + 1 - offset;  % first sample of previous cluster
    x_before(2) = start - 1 - offset; % last sample of previous cluster
    x_after(1)  = start - offset;  % first sample of current (merged) cluster 
    x_after(2)  = stop - offset;  % last sample of current (merged) cluster    
    
end