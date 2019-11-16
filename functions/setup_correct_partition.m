%% Correct partition
function [x, y, x_axis_ticks, nClasses, class_order] = setup_correct_partition(data, classes)
    
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
    
    %% Sort samples    
    class_order = randperm(nClasses);  
    
    x = [];
    y = [];
    x_axis_ticks = zeros(1, nClasses);
    for ix=1:nClasses
        inds = classes == u(class_order(ix));
        nPoints = sum(inds);
        x_temp = data(inds, :);
        % shuffle within cluster
        p = randperm(nPoints);
        x = [x; x_temp(p, :)];
        y = [y; classes(inds)];    
        x_axis_ticks(ix) = find(y == u(class_order(ix)), 1, 'last');
    end   
 
end