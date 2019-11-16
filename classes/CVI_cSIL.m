% PROGRAM DESCRIPTION
% This is a MATLAB implementation of batch and incremental Centroid-based 
% Silhouette Cluster Validity Index.
%
% REFERENCES
% [1] L. E. Brito da Silva, N. M. Melton, and D. C. Wunsch II, "Incremental
% Cluster Validity Indices for Hard Partitions: Extensions  and  Comparative  
% Study," ArXiv  e-prints, Feb 2019, arXiv:1902.06711v1 [cs.LG].
% [2] P. J. Rousseeuw, "Silhouettes: A graphical aid to the interpretation and
% validation of cluster analysis," Journal of Computational and Applied
% Mathematics, vol. 20, pp. 53–65, 1987.
% [3] M. Rawashdeh and A. Ralescu, "Center-wise intra-inter silhouettes," in
% Scalable Uncertainty Management, E. Hüllermeier, S. Link, T. Fober et al.,
% Eds. Berlin, Heidelberg: Springer, 2012, pp. 406–419.
%
% Code written by Leonardo Enzo Brito da Silva
% under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cSIL CVI Class
classdef CVI_cSIL 
    properties (Access = public)       
        dim = [];                   % dimension of the input data  
        nSamples = 0;               % number of samples encountered
        mu = [];                    % geometric mean of the data
        n = [];                     % number of samples belonging to each cluster
        v = [];                     % prototypes (centroids)
        CP = [];                    % compactness of each cluster  
        G = [];                     % vector g of each cluster
        S = [];                     % dissimilarity matrix 
        sil_coefs = [];             % silhouette coefficients
        nClusters = 0;              % number of clusters      
        CriterionValue = [];        % calculated cSIL index
    end 
    methods 
        %% Compute parameters - incremental mode
        function obj = param_inc(obj, sample, label)
            nSamples_new = obj.nSamples + 1;
            if label > obj.nClusters    
                n_new = 1;
                v_new = sample;
                CP_new = sample*sample';
                G_new = sample;  
                if obj.nClusters == 0
                    S_new = 0;
                else
                    S_new = zeros(obj.nClusters + 1, obj.nClusters + 1);
                    S_new(1:obj.nClusters, 1:obj.nClusters) = obj.S;
                    S_row_new = zeros(1, obj.nClusters); 
                    S_col_new = zeros(obj.nClusters, 1);
                    inds = 1:obj.nClusters;                    
                    for cl=inds  
                        % Column "bmu_temp" - D
                        C = CP_new + (obj.v(cl, :)*obj.v(cl, :)') - 2*(obj.v(cl, :)*G_new');
                        S_col_new(cl) = C;                              
                        % Row "bmu_temp" - E
                        C = obj.CP(cl) + obj.n(cl)*(v_new*v_new') - 2*(v_new*obj.G(cl, :)');                       
                        S_row_new(cl) = C/obj.n(cl);
                    end 
                    % Column "ind_minus" - F                          
                    S_col_new(label) = 0;
                    S_row_new(label) = S_col_new(label);
                    
                    S_new(label, :) = S_row_new;
                    S_new(:, label) = S_col_new;
                end   
                % Update parameters
                obj.nClusters = obj.nClusters + 1; 
                obj.n = [obj.n ; n_new];
                obj.v = [obj.v ; v_new]; 
                obj.CP = [obj.CP ; CP_new];     
                obj.G = [obj.G ; G_new]; 
                obj.S = S_new;                
                obj.nSamples = nSamples_new;                   
            else 
                n_new = obj.n(label) + 1; 
                v_new = (1 - 1/n_new).*obj.v(label, :) + (1/n_new).*sample;    
                CP_new = obj.CP(label) + (sample*sample');
                G_new = obj.G(label, :) + sample;
                S_row_new = zeros(1, obj.nClusters); 
                S_col_new = zeros(obj.nClusters, 1);
                inds = 1:obj.nClusters;
                inds(label) = [];
                for cl=inds  
                    % Column "bmu_temp" - D
                    diff_x_v = sample - obj.v(cl, :);
                    C = obj.CP(label) + (diff_x_v*diff_x_v') + obj.n(label)*(obj.v(cl, :)*obj.v(cl, :)') - 2*(obj.v(cl, :)*obj.G(label, :)');
                    S_col_new(cl) = C/n_new;      
                    % Row "bmu_temp" - E
                    C = obj.CP(cl) + obj.n(cl)*(v_new*v_new') - 2*(v_new*obj.G(cl, :)');                                                              
                    S_row_new(cl) = C/obj.n(cl);                                    
                end 
                % Column "ind_minus" - F    
                diff_x_v = sample - v_new;                            
                C = obj.CP(label) + (diff_x_v*diff_x_v') + obj.n(label)*(v_new*v_new') - 2*(v_new*obj.G(label, :)');
                S_col_new(label) = C/n_new;
                S_row_new(label) = S_col_new(label);                
                % Update parameters
                obj.n(label) = n_new;
                obj.v(label, :) = v_new; 
                obj.CP(label) = CP_new;
                obj.G(label, :) = G_new;                
                obj.S(label, :) = S_row_new;
                obj.S(:, label) = S_col_new;
                obj.nSamples = nSamples_new;
            end
        end         
        %% Compute parameters - batch mode
        function obj = param_batch(obj, data, labels)              
            [obj.nSamples, obj.dim] = size(data);
            u = unique(labels);
            obj.nClusters = length(u);   
            obj.n = zeros(obj.nClusters, 1);    
            obj.v = zeros(obj.nClusters, obj.dim);            
            obj.S = zeros(obj.nClusters, obj.nClusters); 
            D = zeros(obj.nClusters, obj.nSamples);
            ind_lbls = cell(1, obj.nClusters);
            for ix=1:obj.nClusters 
                ind_lbls{ix} = labels==u(ix);
                subset = data(ind_lbls{ix}, :); 
                obj.n(ix, 1) = size(subset, 1);
                obj.v(ix, :) = mean(subset, 1);  
                D(ix, 1:obj.nSamples) = (sum((data - ones(obj.nSamples, 1)*obj.v(ix, :)).^2, 2))';        
            end 
            for ix=1:obj.nClusters        
                for jx=1:obj.nClusters            
                    obj.S(ix, jx) = sum(D(ix, ind_lbls{jx}))/obj.n(jx); 
                end
            end              
        end         
        %% Evaluate CVI
        function obj = evaluate(obj)          
                obj.sil_coefs = zeros(obj.nClusters, 1);  
                for ix=1:obj.nClusters                                    
                    ind_other_clusters = true(1, obj.nClusters);
                    ind_other_clusters(ix) = false;
                    a = obj.S(ix, ix);  % same cluster        
                    b = min(obj.S(ix, ind_other_clusters)); % other clusters
                    obj.sil_coefs(ix) = (b - a)/max(a, b);                             
                end 
                obj.CriterionValue = sum(obj.sil_coefs)/obj.nClusters; % cSIL index value 
        end   
    end
end