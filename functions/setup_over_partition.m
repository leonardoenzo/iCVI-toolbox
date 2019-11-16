%% Over-partition
function [x, y, y_over, x_axis_ticks, class_selected, nClasses_over, nClasses, class_order] = setup_over_partition(data, classes)
 
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
    
    %% Complement Coding
    data_CC = [data 1-data];  
    
    %% Compute clusters' sizes 
    R = zeros(1, nClasses);
    nx = zeros(1, nClasses);
    MinPts = 10;
    for ix =1:nClasses
        inds_cl = classes==u(ix);
        nx(ix) = sum(inds_cl);
        if nx(ix) >= MinPts
            W = min(data_CC(inds_cl, :), [], 1);        
            R(ix) = dim - norm(W, 1);
        else
            R(ix) = 0;
        end
    end
    
    %% Select the class (with probability proportional to size)
    class_selected = randsample(nClasses, 1, true, R);

    %% Split the selected class using fuzzy ART in online mode  
      
    % Samples
    inds_selected = classes==class_selected;
    samples_selected = data(inds_selected, :);           
    nPoints_selected = size(samples_selected, 1);      

    % Parameters
    alpha = 1e-3;
    beta = 1;
    rho_large = 1 - R(class_selected)/dim;
    diff_n = Inf;    
    maxEpochs = 1;
    delta = 0.001;
    nTrials = 10;    
    sample_order = []; 
    new_labels = [];
    for tx=1:nTrials
        
        % Complement coded data
        selected_samples_CC = data_CC(inds_selected, :); 
        
        % RNG   
        Prng = randperm(nPoints_selected);
        selected_samples_CC = selected_samples_CC(Prng, :);
        
        % Search loop
        proportion = 1.0;        
        while true
            proportion = proportion + delta;
            rho_small = (1/proportion)*(rho_large + proportion - 1); 
            settings.rho = rho_small;
            settings.alpha = alpha;
            settings.beta = beta; 
            FA = FuzzyART(settings);
            FA.display = false;        
            FA = FA.train(selected_samples_CC, maxEpochs);
            if FA.nCategories>=3
                break;
            elseif FA.nCategories==2  
                n1 = sum(FA.labels==1);
                n2 = sum(FA.labels==2);
                temp_diff_n = abs(n1 - n2);            
                if temp_diff_n <= diff_n  
                    diff_n = temp_diff_n;
                    sample_order = Prng;
                    new_labels = FA.labels + nClasses;
                end
            end         
        end
    end

    % New labels
    classes_over = classes;
    classes_over(inds_selected) = new_labels;    
    um = unique(classes_over);
    nClasses_over = length(um);
    for ix=1:nClasses_over 
        classes_over(classes_over==um(ix)) = ix;
    end
    
    %% Sort samples  
    class_order = randperm(nClasses); 
    if class_order(1) == class_selected  % Avoid the selected cluster to be the first one presented
        pos = datasample(2:nClasses, 1, 'Replace', false);
        var_temp = class_order(pos);
        class_order(pos) = class_order(1);
        class_order(1) = var_temp;
    end
    
    x = [];
    y = [];
    y_over = [];
    x_axis_ticks = zeros(1, nClasses);
    for ix=1:nClasses
        if u(class_order(ix))~=class_selected
            inds2 = classes == u(class_order(ix));
            nPoints2 = sum(inds2);
            x_temp = data(inds2, :);
            % shuffle within cluster
            Prng2 = randperm(nPoints2);
            x = [x; x_temp(Prng2, :)];
            y = [y; classes(inds2)];    
            nINDs = find(inds2==1);
            nINDs = nINDs(Prng2);
            y_over = [y_over; classes_over(nINDs)]; 
            x_axis_ticks(ix) = find(y == u(class_order(ix)), 1, 'last');
        else % samples within the selected cluster were already shuffled            
            x = [x; samples_selected(sample_order,:)];
            y = [y; classes(inds_selected)]; 
            y_over = [y_over; classes_over(inds_selected)]; 
            x_axis_ticks(ix) = find(y == u(class_order(ix)), 1, 'last');
        end
    end  

end