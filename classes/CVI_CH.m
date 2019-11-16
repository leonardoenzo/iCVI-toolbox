% PROGRAM DESCRIPTION
% This is a MATLAB implementation of batch and incremental Calinski-Harabasz (CH)
% Cluster Validity Index.
%
% REFERENCES
% [1] L. E. Brito da Silva, N. M. Melton, and D. C. Wunsch II, "Incremental
% Cluster Validity Indices for Hard Partitions: Extensions  and  Comparative  
% Study," ArXiv  e-prints, Feb 2019, arXiv:1902.06711v1 [cs.LG].
% [2] T. Calinski and J. Harabasz, "A dendrite method for cluster analysis,"
% Communications in Statistics, vol. 3, no. 1, pp. 1–27, 1974.
% [3] M. Moshtaghi, J. C. Bezdek, S. M. Erfani, C. Leckie, and J. Bailey, 
% "Online Cluster Validity Indices for Streaming Data," ArXiv e-prints, 2018, 
% arXiv:1801.02937v1 [stat.ML]. [Online].
% [4] M. Moshtaghi, J. C. Bezdek, S. M. Erfani, C. Leckie, J. Bailey, "Online 
% cluster validity indices for performance monitoring of streaming data clustering," 
% Int. J. Intell. Syst., pp. 1–23, 2018. 
% 
% Code written by Leonardo Enzo Brito da Silva
% under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CH CVI Class
classdef CVI_CH 
    properties (Access = public)       
        dim = [];                   % dimension of the input data  
        nSamples = 0;               % number of samples encountered
        mu = [];                    % geometric mean of the data
        n = [];                     % number of samples belonging to each cluster
        v = [];                     % prototypes (centroids)
        CP = [];                    % compactness of each cluster 
        SEP = [];                   % separation between each cluster and the mean of the data
        G = [];                     % vector g of each cluster
        BGSS = [];                  % between-group sum-of-squares
        WGSS = [];                  % within-group sum-of-squares
        nClusters = 0;              % number of clusters
        CriterionValue = [];        % calculated CH index
    end 
    methods 
        %% Compute parameters - incremental mode
        function obj = param_inc(obj, sample, label)
            nSamples_new = obj.nSamples + 1;
            if isempty(obj.mu)
                mu_new = sample;
            else                
                mu_new = (1 - 1/nSamples_new).*obj.mu + (1/nSamples_new).*sample; 
            end
            if label > obj.nClusters    
                n_new = 1;
                v_new = sample;
                CP_new = 0;
                G_new = zeros(1, obj.dim);  
                % Update parameters
                obj.nClusters = obj.nClusters + 1; 
                obj.n = [obj.n ; n_new];
                obj.v = [obj.v ; v_new]; 
                obj.CP = [obj.CP ; CP_new];      
                obj.G = [obj.G ; G_new]; 
                obj.nSamples = nSamples_new;
                obj.mu = mu_new;   
                for ix=1:obj.nClusters
                    obj.SEP(ix) = obj.n(ix)*sum((obj.v(ix, :) - obj.mu).^2);
                end  
            else 
                n_new = obj.n(label) + 1; 
                v_new = (1 - 1/n_new).*obj.v(label, :) + (1/n_new).*sample;
                Delta_v = obj.v(label, :) - v_new;
                diff_x_v = sample - v_new;     
                CP_new = obj.CP(label) + diff_x_v*(diff_x_v)' + obj.n(label)*Delta_v*(Delta_v)' + 2*Delta_v*obj.G(label, :)';
                G_new = obj.G(label, :) + diff_x_v + obj.n(label).*Delta_v;
                % Update parameters
                obj.n(label) = n_new;
                obj.v(label, :) = v_new; 
                obj.CP(label) = CP_new;
                obj.G(label, :) = G_new;
                obj.nSamples = nSamples_new;
                obj.mu = mu_new;                
                for ix=1:obj.nClusters
                    obj.SEP(ix) = obj.n(ix)*sum((obj.v(ix, :) - obj.mu).^2);
                end  
            end
        end         
        %% Compute parameters - batch mode
        function obj = param_batch(obj, data, labels)  
            [obj.nSamples, obj.dim] = size(data);
            obj.mu = mean(data, 1);
            u = unique(labels);
            obj.nClusters = length(u);   
            obj.n = zeros(obj.nClusters, 1);    
            obj.v = zeros(obj.nClusters, obj.dim);
            obj.CP = zeros(obj.nClusters, 1);
            obj.SEP = zeros(obj.nClusters, 1);
            for ix=1:obj.nClusters
                subset = data(labels==u(ix), :); 
                obj.n(ix) = size(subset, 1);
                obj.v(ix, 1:obj.dim) = mean(subset, 1);
                diff_x_v = subset - ones(obj.n(ix), 1)*obj.v(ix, :);
                obj.CP(ix) = sum(sum(diff_x_v.^2, 2), 1);
                obj.SEP(ix) = obj.n(ix)*sum((obj.v(ix, :) - obj.mu).^2);
            end  
        end         
        %% Evaluate CVI
        function obj = evaluate(obj) 
            obj.WGSS = sum(obj.CP);  
            if obj.nClusters > 1
                obj.BGSS = sum(obj.SEP);                       
                obj.CriterionValue = (obj.BGSS/obj.WGSS)*((obj.nSamples - obj.nClusters)/(obj.nClusters - 1)); % CH index value  
            end
        end   
    end
end