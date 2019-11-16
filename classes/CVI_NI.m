% PROGRAM DESCRIPTION
% This is a MATLAB implementation of batch and incremental Negentropy Increment (NI)
% Cluster Validity Index.
%
% REFERENCES
% [1] L. E. Brito da Silva, N. M. Melton, and D. C. Wunsch II, "Incremental
% Cluster Validity Indices for Hard Partitions: Extensions  and  Comparative  
% Study," ArXiv  e-prints, Feb 2019, arXiv:1902.06711v1 [cs.LG].
% [2] L. F. Lago-Fernandez and F. Corbacho, "Normality-based validation for
% crisp clustering," Pattern Recognition, vol. 43, no. 3, pp. 782 – 795,
% 2010.
% [3] L. F. Lago-Fernandez and F. Corbacho, "Using the Negentropy Increment
% to Determine the Number of Clusters." Berlin, Heidelberg: Springer
% Berlin Heidelberg, 2009, pp. 448–455.
% [4] R. O. Duda, P. E. Hart, and D. G. Stork, Pattern Classification, 2nd ed.
% John Wiley & Sons, 2000.
% 
% Code written by Leonardo Enzo Brito da Silva
% under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NI CVI Class
classdef CVI_NI 
    properties (Access = public)       
        dim = [];                   % dimension of the input data     
        nSamples = 0;               % number of samples encountered
        mu_data = [];               % geometric mean of the data
        Sigma_data = [];            % covariance matrix of the data
        n = [];                     % number of samples belonging to each cluster      
        prob = [];                  % probability of each cluster
        v = [];                     % prototypes (centroids)
        Sigma = [];                 % covariance matrices of clusters
        const = 0;                  % NI term associated solely with the covariance matrix of the data
        pld = [];                   % NI term associated solely with the covariance matrices and probabilities of clusters
        delta_term;                 % small term to avoid numerical errors
        nClusters = 0;              % number of clusters
        CriterionValue = [];        % calculated NI index      
    end 
    methods 
        %% Compute parameters - incremental mode
        function obj = param_inc(obj, sample, label)
            nSamples_new = obj.nSamples + 1;
            if nSamples_new > 1
                mu_data_new = (1 - 1/nSamples_new).*obj.mu_data + (1/nSamples_new).*sample; 
                diff_x_mu = sample - obj.mu_data;
                Sigma_data_new = ((nSamples_new - 2)/(nSamples_new - 1))*(obj.Sigma_data - obj.delta_term) + (1/nSamples_new)*(diff_x_mu'*diff_x_mu);                              
                Sigma_data_new = Sigma_data_new + obj.delta_term; 
            else
                mu_data_new = sample;
                Sigma_data_new = obj.delta_term;
            end             
            if label > obj.nClusters    
                n_new = 1;
                v_new = sample;               
                Sigma_new = obj.delta_term;     
                % Update parameters
                obj.nClusters = obj.nClusters + 1;                
                obj.n = [obj.n ; n_new];
                obj.v = [obj.v ; v_new];                 
                obj.Sigma = cat(3, obj.Sigma, Sigma_new);  
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
                % Update parameters
                obj.n(label) = n_new;
                obj.v(label, :) = v_new;                 
                obj.Sigma(:, :, label) = Sigma_new;            
            end
            obj.prob = obj.n./nSamples_new;
            obj.nSamples = nSamples_new;
            obj.mu_data = mu_data_new;            
            obj.Sigma_data = Sigma_data_new;
        end         
        %% Compute parameters - batch mode
        function obj = param_batch(obj, data, labels)
            [obj.nSamples, obj.dim] = size(data);
            if obj.nSamples > 1
                obj.mu_data = mean(data, 1);
                obj.Sigma_data = (1/(obj.nSamples-1))*((data'*data)-obj.nSamples.*(obj.mu_data'*obj.mu_data)) + obj.delta_term;                
            else
                obj.mu_data = data;
                obj.Sigma_data = obj.delta_term;
            end
            u = unique(labels);
            obj.nClusters = length(u);  
            obj.n = zeros(obj.nClusters, 1);
            obj.prob = zeros(obj.nClusters, 1);
            obj.v = zeros(obj.nClusters, obj.dim);
            obj.Sigma = zeros(obj.dim, obj.dim, obj.nClusters); 
            obj.pld = zeros(obj.nClusters, 1);            
            for ix=1:obj.nClusters                
                subset = data(labels==u(ix), :); 
                obj.n(ix) = size(subset, 1);
                obj.prob(ix) = obj.n(ix)/obj.nSamples;                
                obj.v(ix, 1:obj.dim) = mean(subset, 1);              
                if obj.n(ix) > 1                    
                    obj.Sigma(:, :, ix) = (1/(obj.n(ix) - 1))*((subset'*subset) - obj.n(ix).*(obj.v(ix, :)'*obj.v(ix, :))) + obj.delta_term;                
                else
                    obj.Sigma(:, :, ix) = obj.delta_term;
                end                
            end   
        end         
        %% Evaluate CVI
        function obj = evaluate(obj) 
            for ix=1:obj.nClusters 
                obj.pld(ix) = obj.prob(ix)*log(sqrt(det(obj.Sigma(:, :, ix)))/obj.prob(ix));
            end      
            obj.const = 0.5*log(det(obj.Sigma_data));
            obj.CriterionValue = sum(obj.pld) - obj.const; % NI index value  
        end   
    end
end