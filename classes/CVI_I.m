% PROGRAM DESCRIPTION
% This is a MATLAB implementation of batch and incremental I and Pakhira-Bandyopadhyay-Maulik 
% (PBM) Cluster Validity Indices.
%
% REFERENCES
% [1] L. E. Brito da Silva, N. M. Melton, and D. C. Wunsch II, "Incremental
% Cluster Validity Indices for Hard Partitions: Extensions  and  Comparative  
% Study," ArXiv  e-prints, Feb 2019, arXiv:1902.06711v1 [cs.LG].
% [2] S. Bandyopadhyay and U. Maulik, "Nonparametric genetic clustering:
% comparison of validity indices," IEEE Transactions on Systems, Man,
% and Cybernetics, Part C (Applications and Reviews), vol. 31, no. 1, pp.
% 120–125, Feb 2001.
% [3] M. K. Pakhira, S. Bandyopadhyay, and U. Maulik, "Validity index for
% crisp and fuzzy clusters," Pattern Recognition, vol. 37, no. 3, pp. 487 –
% 501, 2004.
% [4] M. Moshtaghi, J. C. Bezdek, S. M. Erfani, C. Leckie, and J. Bailey, 
% "Online Cluster Validity Indices for Streaming Data," ArXiv e-prints, 2018, 
% arXiv:1801.02937v1 [stat.ML]. [Online].
% [5] M. Moshtaghi, J. C. Bezdek, S. M. Erfani, C. Leckie, J. Bailey, "Online 
% cluster validity indices for performance monitoring of streaming data clustering," 
% Int. J. Intell. Syst., pp. 1–23, 2018. 
% 
% Code written by Leonardo Enzo Brito da Silva
% under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% I and PBM CVI Class
classdef CVI_I 
    properties (Access = public)       
        dim = [];                   % dimension of the input data  
        nSamples = 0;               % number of samples encountered
        mu = [];                    % geometric mean of the data
        n = [];                     % number of samples belonging to each cluster
        v = [];                     % prototypes (centroids)     
        D = [];                     % dissimilarity matrix 
        E = [];                     % compactness of each cluster 
        E0 = [];                    % compactness of the entire data set
        G = [];                     % vector g of each cluster
        G0 = [];                    % vector g of the entire data set    
        p = 2;                      % exponent parameter  (p=2 implies I = PBM)
        nClusters = 0;              % number of clusters
        CriterionValue = [];        % calculated I/PBM index
    end 
    methods 
        %% Compute parameters - incremental mode
        function obj = param_inc(obj, sample, label)
            nSamples_new = obj.nSamples + 1;
            if isempty(obj.mu)
                mu_new = sample;
                E0_new = 0;
                G0_new = zeros(1, obj.dim); 
            else                
                mu_new = (1 - 1/nSamples_new).*obj.mu + (1/nSamples_new).*sample; 
                Delta_v = obj.mu - mu_new;
                diff_x_v = sample - mu_new; 
                E0_new = obj.E0 + diff_x_v*(diff_x_v)' + obj.nSamples*Delta_v*(Delta_v)' + 2*Delta_v*obj.G0';
                G0_new = obj.G0 + diff_x_v + obj.nSamples.*Delta_v; 
            end
            if label > obj.nClusters    
                n_new = 1;
                v_new = sample;
                E_new = 0;                
                G_new = zeros(1, obj.dim); 
                if isempty(obj.v)
                    d_new = [];
                else
                    d_new = sum((obj.v - ones(obj.nClusters, 1)*v_new).^2, 2);
                end              
                % Update parameters
                obj.nClusters = obj.nClusters + 1; 
                obj.n = [obj.n ; n_new];
                obj.v = [obj.v ; v_new]; 
                obj.E = [obj.E ; E_new];                       
                obj.G = [obj.G ; G_new]; 
                obj.D = [obj.D ; d_new'];
                obj.D = [obj.D [d_new ; 0]];                
            else 
                n_new = obj.n(label) + 1; 
                v_new = (1 - 1/n_new).*obj.v(label, :) + (1/n_new).*sample;
                Delta_v = obj.v(label, :) - v_new;
                diff_x_v = sample - v_new;     
                E_new = obj.E(label) + diff_x_v*(diff_x_v)' + obj.n(label)*Delta_v*(Delta_v)' + 2*Delta_v*obj.G(label, :)';
                G_new = obj.G(label, :) + diff_x_v + obj.n(label).*Delta_v;               
                d_new = sum((obj.v - ones(obj.nClusters, 1)*v_new).^2, 2);
                % Update parameters
                obj.n(label) = n_new;
                obj.v(label, :) = v_new; 
                obj.E(label) = E_new;
                obj.G(label, :) = G_new; 
                obj.D(:, label) = d_new;
                obj.D(label, :) = d_new;
                obj.D(label, label) = 0;  
            end
            obj.E0 = E0_new;
            obj.G0 = G0_new;   
            obj.mu = mu_new;
            obj.nSamples = nSamples_new;
        end         
        %% Compute parameters - batch mode
        function obj = param_batch(obj, data, labels) 
            [obj.nSamples, obj.dim] = size(data);
            obj.mu = mean(data, 1);
            u = unique(labels);
            obj.nClusters = length(u);
            obj.n = zeros(obj.nClusters, 1);
            obj.v = zeros(obj.nClusters, obj.dim);            
            obj.E = zeros(obj.nClusters, 1);
            obj.D = zeros(obj.nClusters, obj.nClusters);           
            for ix=1:obj.nClusters
                subset = data(labels==u(ix), :); 
                obj.n(ix) = size(subset, 1);
                obj.v(ix, 1:obj.dim) = mean(subset, 1);                
                diff_x_v = subset - ones(obj.n(ix), 1)*obj.v(ix, :);
                obj.E(ix) = sum(sum(diff_x_v.^2, 2), 1);            
            end
            for ix=1:obj.nClusters-1
                for jx=ix+1:obj.nClusters
                    Delta_v = obj.v(ix, :) - obj.v(jx, :);                    
                    obj.D(ix, jx) = Delta_v*(Delta_v)';  
                end
            end      
            obj.D = obj.D + obj.D';
            diff_x_v = data - ones(obj.nSamples, 1)*obj.mu;
            obj.E0 = sum(sum(diff_x_v.^2, 2), 1);  
        end         
        %% Evaluate CVI
        function obj = evaluate(obj) 
            obj.CriterionValue = ((obj.E0/sum(obj.E))*(max(obj.D(:))/obj.nClusters))^obj.p; % I/PBM index value 
        end   
    end
end