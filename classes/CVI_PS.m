% PROGRAM DESCRIPTION
% This is a MATLAB implementation of batch and incremental Partition Separation (PS)
% Cluster Validity Index.
%
% REFERENCES
% [1] Miin-Shen Yang and Kuo-Lung Wu, "A new validity index for fuzzy clustering," 
% 10th IEEE International Conference on Fuzzy Systems. (Cat. No.01CH37297), Melbourne,
% Victoria, Australia, 2001, pp. 89-92, vol.1.
% [2] E. Lughofer, "Extensions of vector quantization for incremental clustering," Pattern
% Recognit., vol. 41, no. 3, pp. 995–1011, 2008.
% 
% Code written by Leonardo Enzo Brito da Silva
% under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PS CVI Class
classdef CVI_PS 
    properties (Access = public)       
        dim = [];                   % dimension of the input data 
        nSamples = 0;               % number of samples encountered
        n = [];                     % number of samples belonging to each cluster
        v = [];                     % prototypes (centroids)
        D = [];                     % dissimilarity matrix 
        v_bar = [];                 % mean of prototypes
        beta_t = [];                % beta term
        PS_i = [];                  % partition separation term of cluster i
        nClusters = 0;              % number of clusters
        CriterionValue = [];        % calculated PS index    
    end 
    methods 
        %% Compute parameters - incremental mode
        function obj = param_inc(obj, sample, label)
            nSamples_new = obj.nSamples + 1;
            if label > obj.nClusters    
                n_new = 1;
                v_new = sample;
                if obj.nClusters == 0
                    D_new = 0;
                else                
                    D_new = zeros(obj.nClusters + 1, obj.nClusters + 1);
                    D_new(1:obj.nClusters, 1:obj.nClusters) = obj.D;
                    inds = 1:obj.nClusters;
                    d_row_new = zeros(1, obj.nClusters + 1);
                    for jx=inds   
                        d_row_new(jx) = sum((v_new - obj.v(jx, :)).^2, 2);
                    end    
                    D_new(label, :) = d_row_new;
                    D_new(:, label) = d_row_new';  
                end                
                % Update parameters
                obj.nClusters = obj.nClusters + 1; 
                obj.n = [obj.n ; n_new];
                obj.v = [obj.v ; v_new]; 
                obj.D = D_new;
            else 
                n_new = obj.n(label) + 1; 
                v_new = (1 - 1/n_new).*obj.v(label, :) + (1/n_new).*sample;
                inds = 1:obj.nClusters;
                inds(label) = [];
                d_row_new = zeros(1, obj.nClusters);
                for jx=inds   
                    d_row_new(jx) = sum((v_new - obj.v(jx, :)).^2, 2); 
                end                  
                % Update parameters
                obj.n(label) = n_new;
                obj.v(label, :) = v_new;     
                obj.D(label, :) = d_row_new;
                obj.D(:, label) = d_row_new';
            end
            obj.nSamples = nSamples_new;           
        end         
        %% Compute parameters - batch mode
        function obj = param_batch(obj, data, labels)  
            [obj.nSamples, obj.dim] = size(data);
            u = unique(labels);
            obj.nClusters = length(u);   
            obj.n = zeros(obj.nClusters, 1);    
            obj.v = zeros(obj.nClusters, obj.dim);
            obj.D = zeros(obj.nClusters, obj.nClusters);
            obj.PS_i = zeros(obj.nClusters, 1);
            for ix=1:obj.nClusters
                subset = data(labels==u(ix), :); 
                obj.n(ix) = size(subset, 1);
                obj.v(ix, 1:obj.dim) = mean(subset, 1);         
            end   
            for ix=1:obj.nClusters-1
                for jx=ix+1:obj.nClusters
                    Delta_v = obj.v(ix, :) - obj.v(jx, :);                    
                    obj.D(ix, jx) = Delta_v*(Delta_v)';  
                end
            end      
            obj.D = obj.D + obj.D';
        end         
        %% Evaluate CVI
        function obj = evaluate(obj)            
            if obj.nClusters > 1  
                obj.v_bar = mean(obj.v, 1);
                obj.beta_t = 0;
                for ix=1:obj.nClusters
                    Delta_v = obj.v(ix, :) - obj.v_bar;
                    obj.beta_t =  obj.beta_t + Delta_v*(Delta_v)';
                end  
                obj.beta_t = obj.beta_t/obj.nClusters;
                n_max = max(obj.n);
                for ix=1:obj.nClusters
                    d = obj.D(ix, :);
                    d(ix) = []; 
                    obj.PS_i(ix) = (obj.n(ix)/n_max) - exp(-min(d)/obj.beta_t);
                end             
                obj.CriterionValue = sum(obj.PS_i);  % PS index value 
            end
        end   
    end
end