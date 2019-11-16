% PROGRAM DESCRIPTION
% This is a MATLAB implementation of the modified Simplified Fuzzy ARTMAP 
% (SFAM) neural network.
%
% REFERENCES
% [1] L. E. Brito da Silva, N. M. Melton, and D. C. Wunsch II, "Incremental
% Cluster Validity Indices for Hard Partitions: Extensions  and  Comparative  
% Study," ArXiv  e-prints, Feb 2019, arXiv:1902.06711v1 [cs.LG].
% [2] G. A. Carpenter, S. Grossberg, and D. B. Rosen, ''Fuzzy ART: Fast stable
% learning and categorization of analog patterns by an adaptive resonance
% system,'' Neural Networks, vol. 4, no. 6, pp. 759 – 771, 1991.
% [3] G. A. Carpenter, S. Grossberg, N. Markuzon, J. H. Reynolds, and D. B.
% Rosen, ''Fuzzy ARTMAP: A neural network architecture for incremental
% supervised learning of analog multidimensional maps,'' IEEE Transactions
% on Neural Networks, vol. 3, no. 5, pp. 698–713, Sep 1992.
% [4] T. Kasuba, ''Simplified Fuzzy ARTMAP,'' AI Expert, vol. 8, no. 11, pp.
% 18–25, Nov 1993.
%
% Code written by Leonardo Enzo Brito da Silva and Niklas M. Melton 
% Under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modified Simplified Fuzzy ARTMAP
classdef SFuzzyARTMAPMod    
    properties (Access = public)        % default properties' values are set
        ART_A;                          % ART A side
        ART_B;                          % ART B side
        map = []                        % the binary, surjective mapping from B-side to A-side clusters
        labels = [];                    % best matching units (class labels) 
        protos = [];                    % first and second best matching prototype units (A-side)
        dim = []                        % number of variables (columns) in the input data
        nCategories = 0;                % total number of categories
        Epoch = 0;                      % current epoch    
        display = true;                 % displays training progress on the command window (displays intermediate steps)
        W;                              % weight vectors
        n = [];
    end 
    methods        
        % Assign property values from within the class constructor
        function obj = SFuzzyARTMAPMod(settings)
            obj.ART_A = FuzzyART(settings.ART_A); 
            obj.ART_A.display = false;            
        end     
        % Resonance Check 
        function is_res = resonance_check(obj, a, i)
            b = obj.labels(end);
            is_res = true;
            nb = obj.n(b);
            if nb == 2 % force module A's split for every 2nd sample encoded by module B's category
                is_res = false;
            elseif a <= size(obj.map, 1) && any(obj.map(a,[1:b-1 b+1:end]))
                is_res = false;
            end
        end        
        % Train
        function obj = train(obj, x, maxEpochs, y)    
            
            % Display progress on command window
            if obj.display
                fprintf('Starting Training...\n');                        
            end 
            
            % Data Information            
            [nSamples, obj.dim] = size(x);            
            obj.dim = obj.dim/2; 
                      
            % B-side training is independent (supervised: y)              
            if nSamples == 1
                obj.labels = [obj.labels; y];
            else
                obj.labels = y;
                obj.protos = [];
            end            
            obj.nCategories = y;
            if length(obj.n) < y
                obj.n(y) = 1;
            else
                obj.n(y) = obj.n(y) + 1;
            end
            
            % A-side            
            obj.Epoch = 0;
            while(true)
                obj.Epoch = obj.Epoch + 1;
                for i=1:nSamples  % loop over samples
                    % create secondary vigilance check  
                    obj.ART_A.vc2 = []; % prevent matlab's memory issue
                    obj.ART_A.vc2 = @(a) obj.resonance_check(a, i); 
                    obj.ART_A = obj.ART_A.train(x(i, :), 1);
                    obj.protos = [obj.protos; [obj.ART_A.labels, obj.ART_A.labels2]];
                    obj.map(obj.ART_A.labels, y(i)) = 1;
                    % Display Progress on Console Window
                    if obj.display
                       progress = sprintf('\tEpoch: %d \n\tSample ID: %d \n\tCategories: %d \n', obj.Epoch, i, obj.nCategories);
                       fprintf([backspace, progress]);
                       backspace = repmat(sprintf('\b'), 1, length(progress)); 
                    end 
                end 
                % Stopping Conditions
                if stopping_conditions(obj, maxEpochs)
                    break;
                end                                
            end  
            % Display progress on command window
            if obj.display
                fprintf('Done.\n');
            end
        end 
        % Classify
        function [labels_B, labels_A] = classify(obj, x)  
            % B-side
            labels_B = [];
            % A-side 
            [lbls1, lbls2] = obj.ART_A.classify(x);
            labels_A = [lbls1, lbls2];                              
        end       
        % Stopping Criteria
        function stop = stopping_conditions(obj, maxEpochs)
            stop = obj.ART_A.stopping_conditions(maxEpochs, obj.Epoch);                                      
        end    
    end  
end