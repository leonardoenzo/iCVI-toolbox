% PROGRAM DESCRIPTION
% This is a MATLAB implementation of batch and incremental CONN 
% Cluster Validity Index.
%
% REFERENCES
% [1] L. E. Brito da Silva, N. M. Melton, and D. C. Wunsch II, "Incremental
% Cluster Validity Indices for Hard Partitions: Extensions  and  Comparative  
% Study," ArXiv  e-prints, Feb 2019, arXiv:1902.06711v1 [cs.LG].
% [2] K. Tasdemir and E. Merenyi, "A new cluster validity index for prototype
% based clustering algorithms based on inter- and intra-cluster density," in
% Proc. Int. Joint Conf. Neural Netw. (IJCNN), Aug. 2007, pp. 2205–2211.
% [3] K. Tasdemir and E. Merenyi, "A Validity Index for Prototype-Based
% Clustering of Data Sets With Complex Cluster Structures,” IEEE Trans.
% Syst., Man, Cybern. B, vol. 41, no. 4, pp. 1039–1053, Aug. 2011.
%
% Code written by Niklas M. Melton and Leonardo Enzo Brito da Silva
% Under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CONN CVI Class
classdef CVI_CONN 
    properties (Access = public)        
        dim = [];                       % dimension of the input data
        nSamples = 0;                   % number of samples encountered
        label_protos = {};              % prototypes for each cluster
        CADJ = [];                      % CADJ matrix
        CONN = [];                      % CONN matrix
        nClusters = 0;                  % number of clusters
        nPrototypes = 0;                % number of prototypes
        inter_conn = 0;                 % inter connectivity
        intra_conn = 0;                 % intra connectivity
        CriterionValue = [];            % calculated conn_index
        inter_k_cache = [];             % cached average inter conn for each cluster
        inter_kl_cache = [];            % cached inter conn between each cluster
        intra_k_cache = [];             % cahced intra conn for each cluster
        condition = 'CONN';             % experimental switch for interconnectivity (string: 'CONN' or 'CADJ')
        missing_samples = 0;            % track samples with only a single bmu      
    end 
    methods 
        %% Compute parameters - incremental mode
        function obj = param_inc(obj, p, p2, label)
            obj.nSamples = obj.nSamples + 1;                                % increment sample count  
            new_proto_flag = false;                                         % new prototype flag default
            if max(p, p2) > obj.nPrototypes                                 % if new prototype
                new_proto_flag = true;                                      % new prototype flag raised
                obj.nPrototypes = obj.nPrototypes + 1;                      % increment prototype count  
                if p2 > 0                                                   % if more than one prototypes
                    obj.CADJ(p, p2) = 1;                                    % expand CADJ matrix
                    obj.CADJ(p2, p) = 0;
                else
                    obj.CADJ(p, 1) = 0;                                     % initialize CADJ matrix
                    obj.CADJ(1, p) = 0;
                    obj.missing_samples = obj.missing_samples + 1;          % cant track samples in CADJ so must record later
                end
            else
                if p ~= p2 && p2 > 0                                        % if more than one prototype
                    obj.CADJ(p, p2) = obj.CADJ(p, p2) + 1;                  % increment CADJ counter
                end
            end
            if p2 > 0                                                       % if more than one prototype
                if obj.missing_samples > 0                                  % if we have missing samples in the CADJ matrix
                    obj.CADJ(p2, p) = obj.CADJ(p2, p)+obj.missing_samples;  % count samples in CADJ matrix
                    obj.missing_samples = 0;                                % no more missing samples
                end  
                obj.CONN(p, p2) = obj.CADJ(p, p2) + obj.CADJ(p2, p);        % update CONN matrix
                obj.CONN(p2, p) = obj.CONN(p, p2);
            end
            if label > obj.nClusters                                        % if new cluster
                obj.nClusters = obj.nClusters + 1;                          % update cluster count and prototype list                
                obj.label_protos{label} = p;
            elseif new_proto_flag %~ismember(p, obj.label_protos{label})    % if new prototype for this cluster
                obj.label_protos{label} = [obj.label_protos{label}; p];     % remember this cluster-prototype pair
            end  
            obj = obj.calc_intra_conn(label);                               % update cluster intra_conn
            if p2 > 0                                                       % if more than one prototype  
                for Cl=1:obj.nClusters                                     
                    if any(p2==obj.label_protos{Cl}) %ismember(p2, obj.label_protos{Cl})
                        obj = obj.calc_inter_conn(label, Cl);               % update inter_conn for this cluster
                    end
                end
            end            
        end         
        %% Compute parameters - batch mode
        function obj = param_batch(obj, P, P2, labels)            
            
            obj.nSamples = length(labels);                                  % update sample count       
            u_c = unique(labels);
            obj.nClusters = length(u_c);                                    % update cluster count
            u_p = max([P; P2]);                                             % take into account all prototypes, even if not activated
            obj.nPrototypes = u_p;                                          % update prototype count            
            obj.label_protos = cell(obj.nClusters, 1);               
            obj.CADJ = zeros(obj.nPrototypes, obj.nPrototypes);  
            
            for ix=1:obj.nClusters
                obj.label_protos{u_c(ix)} = unique(P(labels==u_c(ix)));
            end
                        
            for i=1:obj.nSamples                                            % for each sample
                if P2(i) == 0 
                    if obj.nPrototypes > P(i)
                        obj.CADJ(P(i), P(i)+1) = obj.CADJ(P(i), P(i)+1) + 1;
                    end
                elseif P(i) ~= P2(i)                                        % increment CADJ matrix
                    obj.CADJ(P(i), P2(i)) = obj.CADJ(P(i), P2(i)) + 1;
                end
            end
            obj.CONN = obj.CADJ + obj.CADJ';                                % build CONN matrix
            obj = obj.calc_intra_conn();                                    % batch calculate intra conn
            obj = obj.calc_inter_conn();                                    % batch calculate inter conn
        end
        %% intra conn for a cluster Ck
        function obj = calc_intra_conn(obj, Ck)  
            if nargin == 2                                                  % if a cluster is given:
                obj.intra_k_cache(Ck) = obj.intra_k(Ck);                    % update cached value for cluster   
            end
            obj.intra_conn = 0;
            for k=1:obj.nClusters
                if nargin == 1                                              % if no cluster is specified
                    obj.intra_k_cache(k) = obj.intra_k(k);                  % update all cached cluster values
                end
                obj.intra_conn = obj.intra_conn + obj.intra_k_cache(k);     % sum cached values
            end
            obj.intra_conn = obj.intra_conn/obj.nClusters;                  % average cached values 
        end
        %% intra conn for a cluster Ck
        function ic = intra_k(obj, Ck)
            protos_k = obj.label_protos{Ck};            
            % Compute numerator
            temp1 = obj.CADJ(protos_k, protos_k);                           % all samples in which 1st & 2nd BMUs belong to Ck
            ic1 = sum(temp1(:)); 
            % Compute denominator                   
            temp2 = obj.CADJ(protos_k, :);                                  % all samples in which 1st BMUs belong to Ck
            ic2 = sum(temp2(:));  
            if ic2 > 0                                                      % return 0 if invalid div
                ic =  ic1/ic2 ;
            else
                ic = 0;
            end
        end            
        %% inter conn between two clusters Ck and Cl
        function obj = calc_inter_conn(obj, Ck, Cl)
            if nargin == 3                                                  % if clusters are specified
                obj = obj.calc_inter_k(Ck, Cl);                             % update Ck, Cl and their connected clusters Cj's (if any)
            else                                                            % ------OTHERWISE------    
                for k=1:obj.nClusters                                       % update all cached values in the inter_conn matrix 
                    obj = obj.calc_inter_k(k);                                                           
                end
            end
            obj.inter_conn = sum(obj.inter_k_cache)/obj.nClusters;          % average inter_conn
        end        
        %% inter Ck for a single cluster Ck or between two clusters Ck,Cl
        function obj = calc_inter_k(obj, Ck, Cl)
            if nargin == 3                                                  % if Ck and Cl clusters are both specified                
                if Ck~=Cl                                                   % only affects clusters Ck and Cl
                    obj.inter_kl_cache(Ck, Cl) = obj.inter_kl(Ck, Cl);
                    obj.inter_kl_cache(Cl, Ck) = obj.inter_kl(Cl, Ck);                    
                    obj.inter_k_cache(Ck) = max(obj.inter_kl_cache(Ck, :));
                    obj.inter_k_cache(Cl) = max(obj.inter_kl_cache(Cl, :));
                else                                                        % only affects inter of Ck
                    for Cj=1:obj.nClusters                                      
                        obj.inter_kl_cache(Ck, Cj) = obj.inter_kl(Ck, Cj);  
                    end
                    obj.inter_k_cache(Ck) = max(obj.inter_kl_cache(Ck, :));
                end
            else                                                            % ------OTHERWISE------     
                for Cl=1:obj.nClusters                                      % update all Ck inter_conn
                    obj.inter_kl_cache(Ck, Cl) = obj.inter_kl(Ck, Cl);
                end
                obj.inter_k_cache(Ck) = max(obj.inter_kl_cache(Ck, :));
            end            
        end        
        %% inter connectivity between two clusters Ck, Cl
        function ic = inter_kl(obj, Ck, Cl)
            protos_k = obj.label_protos{Ck};	% collect protos in Ck
            protos_l = obj.label_protos{Cl};	% collect protos in Cl            
            if Cl ~= Ck  
                % Identify prototypes that belong to cluster borders
                if strcmp(obj.condition, 'CONN')                                    % using CONN(i,j) > 0 as the condition 
                    protos_pkl = protos_k(any(obj.CONN(protos_k, protos_l), 2));    % ----OR----                                                                    % ----OR----
                elseif strcmp(obj.condition, 'CADJ')                                % using CADJ(i,j) > 0 as the condition
                    protos_pkl = protos_k(any(obj.CADJ(protos_k, protos_l), 2));                        
                end                 
                if isempty(protos_pkl)
                    ic = 0;                    
                else
                    % Compute numerator
                    CONN_kl = obj.CONN(protos_pkl, protos_l);
                    conn_ck_cl = sum(CONN_kl(:)); 
                    % Compute denominator 
                    CONN_pkl = obj.CONN(protos_pkl, [protos_k; protos_l]);
                    intersection = obj.CONN(protos_pkl, protos_pkl);                % AND(protos_pkl,[protos_k; protos_l]) = protos_pkl                    
                    conn_pkl = sum(CONN_pkl(:)) - sum(intersection(:))/2;           % subtract what has been counted twice (CONN is symmetric)                  
                    % Compute inter_conn
                    if conn_pkl > 0     % return quotient or zero if invalid div
                        ic = conn_ck_cl / conn_pkl;
                    else
                        ic = 0;
                    end
                end  
            else
                ic = 0;
            end
        end        
        %% Evaluate CVI
        function obj = evaluate(obj)
            obj.CriterionValue = obj.intra_conn*(1-obj.inter_conn);         % Conn_Index value 
        end   
    end
end