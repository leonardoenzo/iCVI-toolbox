%% Under-partition
function [x, y, y_under, x_axis_ticks, classes_merged, nClasses_under, nClasses, class_order] = setup_under_partition(data, classes)
    
    %% Data 
    u = unique(classes);
    nClasses = length(u);
    [nSamples, dim] = size(data);
    
    %% Linear Normalization  
    data = mapminmax(data', 0, 1);
    data = data';     
    for ix=1:dim
        u_values = unique(data(:, ix));
        if length(u_values)==1  % constant 
            data(:, ix) = ones(nSamples, 1); % design choice: set attributes to 1's
        end        
    end      

    %% Compute clusters' simularity (centroids' distance) 
    centroids = zeros(nClasses, dim);
    D = zeros(nClasses, nClasses); 
    for ix=1:nClasses
        centroids(ix, :) = mean(data(classes==u(ix),:));        
    end
    for ix=1:nClasses-1
        for jx=ix+1:nClasses
            delta_v = centroids(ix, :) - centroids(jx, :);                    
            D(ix, jx) = delta_v*(delta_v)';  
        end
    end      
    D = D + D';     
    D = D + diag(Inf(1, nClasses)); 

    %% Select a pair of classes (with probability proportional to their centroid closeness,i.e., inversely proportional to their distance)
    C = nchoosek(1:nClasses, 2);  % all pairwise combinations
    nPairs = size(C, 1);    % total number of pairs
    S = zeros(nPairs, 1);   % similarity vector
    for ix=1:nPairs
        S(ix) = 1/D(C(ix,1), C(ix,2));
    end
    S = S.^3;
    S = S/sum(S);   % normalize: range [0,1]
    pair_selected = randsample(nPairs, 1, true, S);
    classes_merged = C(pair_selected, :);
      
    %% Merge the selected classes
    % New labels
    classes_under = classes;
    classes_under(classes_under==classes_merged(2)) = classes_merged(1);    
    um = unique(classes_under);
    nClasses_under = length(um);
    for ix=1:nClasses_under 
        classes_under(classes_under==um(ix)) = ix;
    end
    
    %% Sort samples    
    class_order = randperm(nClasses); 
    if sum(ismember(class_order(1:2), classes_merged))==2  % Avoid the merged cluster to be the first one presented
        pos = datasample(3:nClasses, 1, 'Replace', false);
        var_temp = class_order(pos);
        class_order(pos) = class_order(2);
        class_order(2) = var_temp;
    end  
    
    x = [];
    y = [];
    y_under = [];
    x_axis_ticks = zeros(1, nClasses);
    for ix=1:nClasses
        inds = classes == u(class_order(ix));
        nPoints = sum(inds);
        x_temp = data(inds, :);
        % shuffle within cluster
        p = randperm(nPoints);
        x = [x; x_temp(p, :)];
        y = [y; classes(inds)];    
        y_under = [y_under; classes_under(inds)]; 
        x_axis_ticks(ix) = find(y == u(class_order(ix)), 1, 'last');
    end 

end