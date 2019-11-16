% PROGRAM DESCRIPTION
% This is a MATLAB implementation of the "Fuzzy ART (FA)" neural network.
%
% REFERENCES
% [1] G. Carpenter, S. Grossberg, and D. Rosen, "Fuzzy ART: Fast 
% stable learning and categorization of analog patterns by an adaptive 
% resonance system," Neural networks, vol. 4, no. 6, pp. 759–771, 1991.
% 
% Code written by Leonardo Enzo Brito da Silva and Niklas M. Melton 
% Under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fuzzy ART Class
classdef FuzzyART    
    properties (Access = public)        % default properties' values
        rho = 0.5;                      % vigilance parameter: [0,1] 
        alpha = 1e-3;                   % choice parameter 
        beta = 1;                       % learning parameter: (0,1] (beta=1: "fast learning")    
        W = [];                         % top-down weights 
        W_old = [];                     % old top-down weight values
        n = [];                         % instance counting (number of samples encoded per category)
        labels = [];                    % best matching units (class labels)
        labels2 = [];                   % second-best matching units
        dim = [];                       % original dimension of data set  
        nCategories = 0;                % total number of categories
        Epoch = 0;                      % current epoch    
        display = true;                 % displays training progress on the command window (displays intermediate steps)
        mismatch_flag;                  % flag that indicates the creation of a new category
        mismatch_flag_2;                % flag that indicates the no second resonant cluster
        T = [];                         % choice function vector
        T_temp = [];                    % temporary choice function vector 
        M = [];                         % match function vector
        M_temp = [];                    % temporary match function vector
        vc2;                            % second vigilance check function
        find_bmu2 = false;              % flag that indicates the second best matching unit is needed        
    end 
    methods        
        % Assign property values from within the class constructor
        function obj = FuzzyART(settings) 
            obj.rho = settings.rho;
            obj.alpha = settings.alpha;
            obj.beta = settings.beta;
        end         
        % Train
        function obj = train(obj, x, maxEpochs)              
            % Display progress on command window
            if obj.display
                fprintf('Starting Training...\n');                
                backspace = '';          
            end              
            % Data Information            
            [nSamples, obj.dim] = size(x);
            obj.dim = obj.dim/2;
            obj.labels = zeros(nSamples, 1);                        
            % Initialization 
            if isempty(obj.W)             
                obj.W = ones(1, 2*obj.dim);    
                obj.n = 0; 
                obj.nCategories = 1;                 
            end              
            obj.W_old = obj.W;  
            % Learning            
            obj.Epoch = 0;
            while(true)
                obj.Epoch = obj.Epoch + 1;
                for i=1:nSamples  % loop over samples 
                    if or(isempty(obj.T), isempty(obj.M))       % Check for already computed activation/match values
                        obj = activation_match(obj, x(i,:));	% Compute Activation/Match Functions
                    end 
                    [~, index] = sort(obj.T, 'descend');  % Sort activation function values in descending order                    
                    obj.mismatch_flag = true;  % mismatch flag 
                    obj.mismatch_flag_2 = true;  % mismatch flag for second best matching unit 
                    for j=1:obj.nCategories  % loop over categories  
                        bmu = index(j);  % Best Matching Unit 
                        if obj.M(bmu) >= obj.rho*obj.dim % Vigilance Check - Pass 
                            if ~isa(obj.vc2,'function_handle') || obj.vc2(bmu) % secondary vigilance check for ARTMAP
                                if obj.mismatch_flag
                                    obj = learn(obj, x(i,:), bmu);  % learning
                                    obj.labels(i) = bmu;  % update sample labels
                                    obj.mismatch_flag = false;  % mismatch flag 
                                    if ~obj.find_bmu2
                                        obj.mismatch_flag_2 = false; % mismatch flag 2
                                        break;
                                    end
                                else
                                    obj.labels2(i) = bmu;
                                    break; 
                                end
                            end
                        end                               
                    end  
                    if obj.mismatch_flag  % If there was no resonance at all then create new category
                        obj.nCategories = obj.nCategories + 1;   % increment number of categories
                        obj.W(obj.nCategories, :) = x(i, :);     % fast commit 
                        obj.n(obj.nCategories, 1) = 1;           % increment instance counter
                        obj.labels(i) = obj.nCategories;         % update sample labels
                        if obj.find_bmu2
                            obj.labels2(i) = index(1);           % update second label (2nd BMU not tested for resonance)
                            obj.mismatch_flag_2 = false;         % mismatch flag 2
                        end
                        obj.T_temp(end+1) = obj.dim/(obj.alpha + obj.dim);  % store T values in a temporary variable
                        obj.M_temp(end+1) = obj.dim;  % store M values in a temporary variable                         
                    end
                    if obj.find_bmu2 && obj.mismatch_flag_2
                        if obj.nCategories > 1
                            if obj.labels(i) == index(1)
                                obj.labels2(i) = index(2);
                            else  
                                obj.labels2(i) = index(1);
                            end
                        else
                            obj.labels2(i) = 0;
                        end
                    end     
                    obj.T = [];  % empty activation vector
                    obj.M = [];  % empty match vector                    
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
                obj.W_old = obj.W;  
            end 
            % Display progress on command window
            if obj.display
                fprintf('Done.\n');
            end            
        end  
        % Individual activation function values (i.e., per category)
        function Tc = activation(obj, x, index)            
            % Data Information            
            nSamples = size(x, 1);  
            % Allocate memory
            Tc = zeros(1, nSamples);
            % Compute activation functions 
            for ix=1:nSamples  % loop over samples                 
                numerator = norm(min(x(ix, :), obj.W(index, :)), 1);
                Tc(1, ix) = numerator/(obj.alpha + norm(obj.W(index, :), 1));  
            end          
        end          
        % Classify
        function [lbls1, lbls2, Tc] = classify(obj, x)                         
            % Data Information            
            nSamples = size(x, 1);            
            lbls1 = zeros(nSamples, 1);                        
            lbls2 = zeros(nSamples, 1);   
            % Classify
            for ix=1:nSamples  % loop over samples 
                % Compute Activation/Match Functions 
                Tc = zeros(obj.nCategories, 1);     
                Mc = zeros(obj.nCategories, 1); 
                for jx=1:obj.nCategories 
                    numerator = norm(min(x(ix, :), obj.W(jx, :)), 1);
                    Tc(jx, 1) = numerator/(obj.alpha + norm(obj.W(jx, :), 1));
                    Mc(jx, 1) = numerator;
                end	 
                [~, index] = sort(Tc, 'descend');	% Sort activation function values in descending order                    
                lbls1(ix) = index(1);               % update sample labels (1st BMU)
                lbls2(ix) = index(2);               % update sample labels (2nd BMU)      
            end
        end                
        % Activation/Match Functions
        function obj = activation_match(obj, x)              
            obj.T = zeros(obj.nCategories, 1);     
            obj.M = zeros(obj.nCategories, 1); 
            for j=1:obj.nCategories
                numerator = norm(min(x, obj.W(j, :)), 1);
                obj.T(j, 1) = numerator/(obj.alpha + norm(obj.W(j, :), 1));
                obj.M(j, 1) = numerator;
            end
            obj.T_temp = obj.T;  % store T values in a temporary variable
            obj.M_temp = obj.M;  % store M values in a temporary variable
        end  
        % Learning
        function obj = learn(obj, x, index)
            obj.W(index, :) = obj.beta*(min(x, obj.W(index, :))) + (1-obj.beta)*obj.W(index, :); 
            obj.n(index, 1) = obj.n(index, 1) + 1;
        end      
        % Stopping Criteria
        function stop = stopping_conditions(obj, maxEpochs, epochs)
            if nargin < 3
                epochs = obj.Epoch;
            end
            stop = false; 
            if isequal(obj.W, obj.W_old)
                stop = true;                                         
            elseif epochs >= maxEpochs
                stop = true;
            end 
        end    
    end  
end