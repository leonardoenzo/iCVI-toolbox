%% This is an example of the iCVI toolbox usage.
%
% PROGRAM DESCRIPTION
% This program exemplifies the usage of the iCVI toolbox provided. The
% latter includes code for batch and incremental versions of:
% 	1.  Calinski-Harabasz
% 	2.  I index
% 		Additional Parameter: p>=1 (default = 2)
% 		Note: The default reduces the I index to Pakhira-Bandyopadhyay-Maulik
% 	3.  WB index
% 	4.  Silhouette
% 	5.  Xie-Beni
% 	6.  Davies-Bouldin
% 	7.  Generalized Dunn's Index 43
% 	8.  Generalized Dunn's Index 53
% 	9.  Partition Separation
% 	10. Negentropy Increment
% 		Additional Parameter: epsilon (default = 12)
% 	11. Representative Cross Information Potential
% 		Additional Parameter: epsilon (default = 12)
% 	12. Representative Cross Entropy
% 		Additional Parameter: epsilon (default = 12)
% 	13. Conn_Index
% 		Additional Parameters: condition 'CONN' or 'CADJ' (default = 'CONN')
%                              ART-A vigilance parameter rho_A
% 		Note: In this example rho_A = 0.9
% 
% REFERENCES
% [1] L. E. Brito da Silva, N. M. Melton, and D. C. Wunsch II, "Incremental
% Cluster Validity Indices for Hard Partitions: Extensions  and  Comparative  
% Study," ArXiv  e-prints, Feb 2019, arXiv:1902.06711v1 [cs.LG].
%
% Code written by Leonardo Enzo Brito da Silva
% under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean up
clear variables; close all; fclose all; echo off; clc;

%% Add path
addpath('classes', 'functions', 'data', 'inputs');

%% Inputs
fprintf('************************************************************************************\n')
fprintf('************************************************************************************\n')
fprintf('Cluster Validity Indices (CVIs) available:\n')
fprintf('\t1.  Calinski-Harabasz: for this CVI type CH\n')
fprintf('\t2.  I index: for this CVI type I\n')
fprintf('\t3.  WB index: for this CVI type WB\n')
fprintf('\t4.  Silhouette: for this CVI type SIL\n')
fprintf('\t5.  Xie-Beni: for this CVI type XB\n')
fprintf('\t6.  Davies-Bouldin: for this CVI type DB\n')
fprintf('\t7.  Generalized Dunn''s Index 43: for this CVI type GD43\n')
fprintf('\t8.  Generalized Dunn''s Index 53: for this CVI type GD53\n')
fprintf('\t9.  Partition Separation: for this CVI type PS\n')
fprintf('\t10. Negentropy Increment: for this CVI type NI\n')
fprintf('\t11. Representative Cross Information Potential: for this CVI type rCIP\n')
fprintf('\t12. Representative Cross Entropy: for this CVI type rH\n')
fprintf('\t13. Conn_Index: for this CVI type CONN\n')
set_icvis_names();
opts_CVI = get_icvis_names(); 
CVI = input('Please select an iCVI to continue...\n', 's');
ind_CVI = ismember(opts_CVI, CVI);
if any(ind_CVI)
    ind2 = find(ind_CVI);
    switch ind2
        case 1
            CVI_name = 'Calinski-Harabasz';
        case 2
            CVI_name = 'Silhouette';
        case 3        
            CVI_name = 'Pakhira-Bandyopadhyay-Maulik';
        case 4
            CVI_name = 'WB';
        case 5
            CVI_name = 'Xie-Beni';
        case 6
            CVI_name = 'Davies-Bouldin';
        case 7
            CVI_name = 'Generalized Dunn''s Index 43';
        case 8
            CVI_name = 'Generalized Dunn''s Index 53';
        case 9  
            CVI_name = 'Partition Separation';
        case 10
            CVI_name = 'Conn_Index';
        case 11
            CVI_name = 'Negentropy Increment';
        case 12        
            CVI_name = 'Representative Cross Information Potential';
        case 13        
            CVI_name = 'Representative Cross Entropy';
    end
    fprintf('\n\tIncremental CVI selected: %s\n\n', CVI_name);
    
else
    error('\nThe input provided was: %s. \nThe input (CVI) must be one of the following:\nCH, SIL, I, WB, XB, DB, GD43, GD53, PS, CONN, NI, rCIP or rH.', CVI)
end
fprintf('Experiments available:\n')
fprintf('\t1.  Correct-partitioning: for this experiment type cp\n')
fprintf('\t2.  Under-partitioning: for this experiment type up\n')
fprintf('\t3.  Over-partitioning: for this experiment type op\n')
input_experiment = input('Please select an experiment type to continue...\n', 's');
opts_exp = {'cp', 'up', 'op'}; 
ind_exp = ismember(opts_exp, input_experiment);
if any(ind_exp)
    ind2 = find(ind_exp);
    switch ind2
        case 1
        	experiment = 'correct-partition'; 
        case 2
        	experiment = 'under-partition'; 
        case 3
        	experiment = 'over-partition'; 
    end    
    fprintf('\n\tExperiment selected: %s\n\n', experiment);
else
    error('\nThe input provided was: %s. \nThe input (experiment) must be one of the following:\ncp, up or op.', input_experiment)
end
fprintf('************************************************************************************\n')
fprintf('************************************************************************************\n')

%% Load Data
fprintf('Data set selected: D4\n');
fprintf('Loading data...\n');
load D4.mat 
[nSamples, dim] = size(data);  
fprintf('Done.\n');
fprintf('************************************************************************************\n')
fprintf('************************************************************************************\n')

%% CVI  
fprintf('Setting up CVI...\n');
SFAM = [];
switch CVI
    case 'CH'
        valind_inc = CVI_CH(); 
    case 'SIL'
        valind_inc = CVI_cSIL(); 
    case 'I'        
        valind_inc = CVI_I();
        valind_inc.p = 2;  
    case 'WB'
        valind_inc = CVI_WB();
        CVI_name = 'WB';
    case 'XB'
        valind_inc = CVI_XB();
    case 'DB'
        valind_inc = CVI_DB();
     case 'GD43'
        valind_inc = CVI_GD43();
    case 'GD53'
        valind_inc = CVI_GD53();
    case 'NI'  
        valind_inc = CVI_NI();   
    case 'rCIP'        
        valind_inc = CVI_rCIP();
    case 'rH'        
        valind_inc = CVI_rH();   
    case 'PS'
        valind_inc = CVI_PS();  
    case 'CONN'
        valind_inc = CVI_CONN();
        valind_inc.condition = 'CONN';      
        % Modified Simplified Fuzzy ARTMAP
        settings = struct();    
        settings.ART_A = struct();
        settings.ART_A.rho = 0.9;
        settings.ART_A.alpha = 1e-3;
        settings.ART_A.beta = 1; 
        SFAM = SFuzzyARTMAPMod(settings); 
        SFAM.display = false; 
        SFAM.ART_A.find_bmu2 = true;
end
valind_inc.dim = dim;
if strcmp(CVI, 'rH') || strcmp(CVI, 'rCIP') || strcmp(CVI, 'NI')
    epsilon = 12;
    delta = 10^(-epsilon/dim);
    valind_inc.delta_term = delta.*eye(dim, dim);
end
fprintf('Done.\n');
fprintf('************************************************************************************\n')
fprintf('************************************************************************************\n')
  
%% Experimental example
fprintf('Starting experiment...\n');
rng(0, 'twister'); % Set seed for reproducibility 
switch experiment
    case 'correct-partition'
        fprintf('\t\tSetting up correct partition...\n');
        [x, y, x_axis_ticks, nClasses, class_order] = setup_correct_partition(data, classes);
    case 'under-partition' 
        fprintf('\t\tSetting up under-partition...\n');
        [x, y, y_under, x_axis_ticks, classes_merged, nClasses_under, nClasses, class_order] = setup_under_partition(data, classes);
    case 'over-partition'
        fprintf('\t\tSetting up over-partition...\n');
        [x, y, y_over, x_axis_ticks, class_selected, nClasses_over, nClasses, class_order] = setup_over_partition(data, classes);
end
fprintf('\t\tDone.\n');

% Allocate Variables
CriterionValue_inc = zeros(size(x, 1), 1);    
nClusters_inc = zeros(size(x, 1), 1);    
labels = zeros(size(x, 1), 1);
protos = zeros(size(x, 1), 2);        
switch experiment
    case 'correct-partition'
        labels_tracking = zeros(1, nClasses);
        y_hat = y;  
    case 'under-partition' 
        labels_tracking = zeros(1, nClasses_under);
        y_hat = y_under;
    case 'over-partition'
        labels_tracking = zeros(1, nClasses_over);
        y_hat = y_over;
end

% Experiment  
maxEpochs = 1;    
counter = 1;
fprintf('\t\tStreaming samples...\n');
for ix = 1:nSamples
    % Update iCVI parameters (Clustering Algorithm Agnostic)    
    if ix==1
       labels(ix) = 1; 
       labels_tracking(counter) = y_hat(ix);            
    else
        if ~ismember(y_hat(ix), labels_tracking) % new cluster
            labels(ix) = sum(labels_tracking>0) + 1;
            counter = counter + 1;
            labels_tracking(counter) = y_hat(ix);
        else  
            labels(ix) = find(labels_tracking == y_hat(ix)) ;
        end
    end

    % Update iCVI parameters                 
    if ~strcmp(CVI, 'CONN')                    
        valind_inc = valind_inc.param_inc(x(ix, :), labels(ix));
    else % ART
        xcc = [x(ix, :) 1 - x(ix, :)];  % Complemenent-code 
        SFAM = SFAM.train(xcc, maxEpochs, labels(ix));                     
        protos(ix, :) = [SFAM.ART_A.labels, SFAM.ART_A.labels2];
        valind_inc = valind_inc.param_inc(protos(ix, 1), protos(ix, 2), labels(ix));
    end

    % Track variables 01                
    nClusters_inc(ix) = valind_inc.nClusters;

    % Compute iCVI if there are at least 2 clusters
    if nClusters_inc(ix) > 1
        valind_inc = valind_inc.evaluate(); 
    else
        valind_inc.CriterionValue = 0;
    end    

    % Track variables 02
    CriterionValue_inc(ix) = valind_inc.CriterionValue;

end 
fprintf('\t\tDone.\n');
    
%% Plot behaviors/partitions
fprintf('\t\tPlotting outputs...\n');
switch experiment
    case 'correct-partition'
        plot_behavior_correct_partition(x, labels, CVI, CriterionValue_inc, nClusters_inc, valind_inc, SFAM, x_axis_ticks, CVI_name)
    case 'under-partition' 
        plot_behavior_under_partition(x, labels, CVI, CriterionValue_inc, nClusters_inc, valind_inc, SFAM, x_axis_ticks, CVI_name, class_order, classes_merged, y)
    case 'over-partition'
        plot_behavior_over_partition(x, labels, CVI, CriterionValue_inc, nClusters_inc, valind_inc, SFAM, x_axis_ticks, CVI_name, class_selected, y)
end     
fprintf('\t\tDone.\n');
fprintf('Done.\n');
fprintf('************************************************************************************\n')
fprintf('************************************************************************************\n')