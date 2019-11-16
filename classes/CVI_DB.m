% PROGRAM DESCRIPTION
% This is a MATLAB implementation of batch and incremental Davies-Bouldin (DB) 
% Cluster Validity Index.
%
% REFERENCES
% [1] D. L. Davies and D. W. Bouldin, "A cluster separation measure," 
% IEEE Transaction on Pattern Analysis and Machine Intelligence, vol. 1,
% no. 2, pp. 224–227, Feb. 1979.
% [2] M. Moshtaghi, J. C. Bezdek, S. M. Erfani, C. Leckie, and J. Bailey, 
% "Online Cluster Validity Indices for Streaming Data," ArXiv e-prints, 2018, 
% arXiv:1801.02937v1 [stat.ML]. [Online].
% [3] M. Moshtaghi, J. C. Bezdek, S. M. Erfani, C. Leckie, J. Bailey, "Online 
% cluster validity indices for performance monitoring of streaming data clustering," 
% Int. J. Intell. Syst., pp. 1–23, 2018.
% 
% Code written by Leonardo Enzo Brito da Silva
% under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DB CVI Class
classdef CVI_DB 
    properties (Access = public)       
        dim = [];                   % dimension of the input data      
        nSamples = 0;               % number of samples encountered
        mu_data = [];               % geometric mean of the data
        n = [];                     % number of samples belonging to each cluster
        v = [];                     % prototypes (centroids)
        CP = [];                    % compactness of each cluster 
        S = [];                     % a measure of compactness/dispersion of each cluster
        R = [];                     % similarity of each cluster to the most similar (and different) cluster  
        G = [];                     % vector g of each cluster
        D = [];                     % dissimilarity matrix 
        nClusters = 0;              % number of clusters
        CriterionValue = [];        % calculated DB index
    end 
    methods 
        %% Compute parameters - incremental mode
        function obj = param_inc(obj, sample, label)
            nSamples_new = obj.nSamples + 1;
            if isempty(obj.mu_data)
                mu_data_new = sample;
            else                
                mu_data_new = (1 - 1/nSamples_new).*obj.mu_data + (1/nSamples_new).*sample; 
            end
            if label > obj.nClusters    
                n_new = 1;
                v_new = sample;
                CP_new = 0;                
                G_new = zeros(1, obj.dim);  
                S_new = 0; 
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
                obj.CP = [obj.CP ; CP_new];    
                obj.G = [obj.G ; G_new];   
                obj.D = D_new;
                obj.S = [obj.S ; S_new];
            else 
                n_new = obj.n(label) + 1; 
                v_new = (1 - 1/n_new).*obj.v(label, :) + (1/n_new).*sample;
                Delta_v = obj.v(label, :) - v_new;
                diff_x_v = sample - v_new;     
                CP_new = obj.CP(label) + diff_x_v*(diff_x_v)' + obj.n(label)*Delta_v*(Delta_v)' + 2*Delta_v*obj.G(label, :)';
                G_new = obj.G(label, :) + diff_x_v + obj.n(label).*Delta_v;
                S_new = CP_new/n_new;
                inds = 1:obj.nClusters;
                inds(label) = [];
                d_row_new = zeros(1, obj.nClusters);
                for jx=inds   
                    d_row_new(jx) = sum((v_new - obj.v(jx, :)).^2, 2); 
                end                  
                % Update parameters
                obj.n(label) = n_new;
                obj.v(label, :) = v_new; 
                obj.CP(label) = CP_new;
                obj.G(label, :) = G_new; 
                obj.S(label) = S_new;
                obj.D(label, :) = d_row_new;
                obj.D(:, label) = d_row_new';
            end
            obj.nSamples = nSamples_new;
            obj.mu_data = mu_data_new;             
        end         
        %% Compute parameters - batch mode
        function obj = param_batch(obj, data, labels)  
            [obj.nSamples, obj.dim] = size(data);
            obj.mu_data = mean(data, 1);
            u = unique(labels);
            obj.nClusters = length(u);   
            obj.n = zeros(obj.nClusters, 1);    
            obj.v = zeros(obj.nClusters, obj.dim);
            obj.CP = zeros(obj.nClusters, 1); 
            obj.D = zeros(obj.nClusters, obj.nClusters);            
            for ix=1:obj.nClusters
                subset = data(labels==u(ix), :); 
                obj.n(ix) = size(subset, 1);
                obj.v(ix, 1:obj.dim) = mean(subset, 1);
                diff_x_v = subset - ones(obj.n(ix), 1)*obj.v(ix, :);
                obj.CP(ix) = sum(sum(diff_x_v.^2, 2), 1); 
                obj.S(ix) = obj.CP(ix)/obj.n(ix);
            end 
            for ix=1:(obj.nClusters-1)
                for jx=ix+1:obj.nClusters
                    obj.D(ix, jx) = sum((obj.v(ix, :) - obj.v(jx, :)).^2, 2);                    
                end
            end
            obj.D = obj.D + obj.D';          
        end         
        %% Evaluate CVI
        function obj = evaluate(obj) 
            obj.R = zeros(obj.nClusters, obj.nClusters); 
            for ix=1:(obj.nClusters-1)
                for jx=ix+1:obj.nClusters                         
                    obj.R(ix, jx) = (obj.S(ix) + obj.S(jx))/obj.D(ix, jx); 
                end
            end        
            obj.R = obj.R + obj.R';
            obj.CriterionValue = sum(max(obj.R, [], 1), 2)/obj.nClusters;
        end   
    end
end