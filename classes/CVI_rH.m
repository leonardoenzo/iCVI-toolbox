% PROGRAM DESCRIPTION
% This is a MATLAB implementation of batch and incremental (Renyi's) 
% representative Cross Entropy (rH) Cluster Evaluation Function (CEF).
%
% REFERENCES
% [1] L. E. Brito da Silva, N. M. Melton, and D. C. Wunsch II, "Incremental
% Cluster Validity Indices for Hard Partitions: Extensions  and  Comparative  
% Study," ArXiv  e-prints, Feb 2019, arXiv:1902.06711v1 [cs.LG].
% [2] E. Gokcay and J. C. Principe, "A new clustering evaluation function
% using Renyi’s information potential," in Proc. Int. Conf. Acoust., Speech,
% Signal Process. (ICASSP), vol. 6. Jun. 2000, pp. 3490–3493.
% [3] E. Gokcay and J. C. Principe, "Information theoretic clustering," IEEE
% Trans. Pattern Anal. Mach. Intell., vol. 24, no. 2, pp. 158–171, Feb. 2002.
% [4] D. Araújo, A. D. Neto, and A. Martins, "Representative cross information
% potential clustering," Pattern Recognit. Lett., vol. 34, no. 16,
% pp. 2181–2191, Dec. 2013.
% [5] D. Araújo, A. D. Neto, and A. Martins, "Information-theoretic clustering:
% A representative and evolutionary approach," Expert Syst. Appl.,
% vol. 40, no. 10, pp. 4190–4205, Aug. 2013.
% [6] R. O. Duda, P. E. Hart, and D. G. Stork, Pattern Classification, 2nd ed.
% John Wiley & Sons, 2000.
% 
% Code written by Leonardo Enzo Brito da Silva
% under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% rH/CEF CVI Class
classdef CVI_rH 
    properties (Access = public)       
        dim = [];                   % dimension of the input data
        nSamples = 0;               % number of samples encountered
        n = [];                     % number of samples belonging to each cluster
        v = [];                     % prototypes (centroids)
        Sigma = [];                 % covariance matrices of clusters
        D = [];                     % dissimilarity matrix
        delta_term;                 % small term to avoid numerical errors
        nClusters = 0;              % number of clusters
        CriterionValue = [];        % calculated rCIP index
    end 
    methods 
        %% Compute parameters - incremental mode
        function obj = param_inc(obj, sample, label)
            nSamples_new = obj.nSamples + 1;            
            if label > obj.nClusters    
                n_new = 1;
                v_new = sample;
                Sigma_new = obj.delta_term;    
                if obj.nClusters == 0
                    D_new = 0;
                else                
                    D_new = zeros(obj.nClusters + 1, obj.nClusters + 1);
                    D_new(1:obj.nClusters, 1:obj.nClusters) = obj.D;
                    inds = 1:obj.nClusters;
                    d_row_new = zeros(1, obj.nClusters + 1);
                    for jx=inds
                        diff_m = v_new - obj.v(jx, :);
                        Sigma_q = Sigma_new + obj.Sigma(:, :, jx);
                        d_row_new(jx) = 0.5*(diff_m*inv(Sigma_q)*diff_m' + log(det(Sigma_q)) + obj.dim*log(2*pi));
                    end    
                    D_new(label, :) = d_row_new;
                    D_new(:, label) = d_row_new';  
                end
                % Update parameters
                obj.nClusters = obj.nClusters + 1;                
                obj.n = [obj.n ; n_new];
                obj.v = [obj.v ; v_new]; 
                obj.Sigma = cat(3, obj.Sigma, Sigma_new);  
                obj.D = D_new;
            else 
                n_new = obj.n(label) + 1; 
                v_new = (1 - 1/n_new).*obj.v(label, :) + (1/n_new).*sample;
                diff_x_v = sample - obj.v(label, :);
                if n_new > 1
                    Sigma_new = ((n_new - 2)/(n_new - 1))*(obj.Sigma(:, :, label) - obj.delta_term) + (1/n_new)*(diff_x_v'*diff_x_v);
                    Sigma_new = Sigma_new + obj.delta_term; 
                else
                    Sigma_new = obj.delta_term;
                end      
                inds = 1:obj.nClusters;
                inds(label) = [];
                d_row_new = zeros(1, obj.nClusters);
                for jx=inds
                    diff_m = v_new - obj.v(jx, :);
                    Sigma_q = Sigma_new + obj.Sigma(:, :, jx);
                    d_row_new(jx) = 0.5*(diff_m*inv(Sigma_q)*diff_m' + log(det(Sigma_q)) + obj.dim*log(2*pi)); 
                end  
                % Update parameters
                obj.n(label) = n_new;
                obj.v(label, :) = v_new; 
                obj.Sigma(:, :, label) = Sigma_new; 
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
            obj.Sigma = zeros(obj.dim, obj.dim, obj.nClusters); 
            for ix=1:obj.nClusters                
                subset = data(labels==u(ix), :); 
                obj.n(ix) = size(subset, 1);                              
                obj.v(ix, 1:obj.dim) = mean(subset, 1);                
                if obj.n(ix) > 1                    
                    obj.Sigma(:, :, ix) = (1/(obj.n(ix) - 1))*((subset'*subset) - obj.n(ix).*obj.v(ix, :)'*obj.v(ix, :)) + obj.delta_term;
                else
                    obj.Sigma(:, :, ix) = obj.delta_term;
                end                
            end         
            obj.D = zeros(obj.nClusters, obj.nClusters);
            for ix=1:(obj.nClusters-1)
                for jx=ix+1:obj.nClusters
                    diff_m = obj.v(ix, :) - obj.v(jx, :);
                    Sigma_q = obj.Sigma(:, :, ix) + obj.Sigma(:, :, jx);
                    obj.D(ix, jx) = 0.5*(diff_m*inv(Sigma_q)*diff_m' + log(det(Sigma_q)) + obj.dim*log(2*pi));                    
                end
            end            
            obj.D = obj.D + obj.D';
        end         
        %% Evaluate CVI
        function obj = evaluate(obj) 
            mask = ones(size(obj.D));
            mask = logical(triu(mask, 1));
            values = obj.D(mask);
            obj.CriterionValue = sum(values); % rH/CEF index value
        end   
    end
end