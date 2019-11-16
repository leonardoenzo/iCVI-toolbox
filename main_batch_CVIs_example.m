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
% 		Note: In this example rho_A = 0.0
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
    fprintf('\n\tBatch CVI selected: %s\n\n', CVI_name);
else
    error('\nThe input provided was: %s. \nThe input (CVI) must be one of the following:\nCH, SIL, I, WB, XB, DB, GD43, GD53, PS, CONN, NI, rCIP or rH.', CVI)
end

%% Load Data
fprintf('Data set selected: D4\n');
fprintf('Loading and processing data...\n');
load D4.mat 
% Info
[nSamples, dim] = size(data);  
u = unique(classes);
nClasses = length(u);    
% Linear Normalization  
data = mapminmax(data', 0, 1);
data = data';     
for ix=1:dim
    u_values = unique(data(:, ix));
    if length(u_values)==1  % constant 
        data(:, ix) = ones(nSamples, 1); % design choice: set attributes to 1's
    end        
end   
fprintf('Done.\n');
fprintf('************************************************************************************\n')
fprintf('************************************************************************************\n')

%% CVI  
fprintf('Setting up CVI...\n');
ART = [];
switch CVI
    case 'CH'
        valind_batch = CVI_CH(); 
    case 'SIL'
        valind_batch = CVI_cSIL(); 
    case 'I'        
        valind_batch = CVI_I();
        valind_batch.p = 2;  
    case 'WB'
        valind_batch = CVI_WB();
    case 'XB'
        valind_batch = CVI_XB();
    case 'DB'
        valind_batch = CVI_DB();
    case 'GD43'
        valind_batch = CVI_GD43();
    case 'GD53'
        valind_batch = CVI_GD53();
    case 'NI'  
        valind_batch = CVI_NI();   
    case 'rCIP'        
        valind_batch = CVI_rCIP();
    case 'rH'        
        valind_batch = CVI_rH();   
    case 'PS'
        valind_batch = CVI_PS();  
    case 'CONN'
        valind_batch = CVI_CONN();
        valind_batch.condition = 'CONN';    
        % Modified Simplified Fuzzy ARTMAP
        settings = struct();    
        settings.ART_A = struct();
        settings.ART_A.rho = 0.0;
        settings.ART_A.alpha = 1e-3;
        settings.ART_A.beta = 1; 
        SFAM = SFuzzyARTMAPMod(settings); 
        SFAM.display = false; 
        SFAM.ART_A.find_bmu2 = true;
end
valind_batch.dim = dim;
if strcmp(CVI, 'rH') || strcmp(CVI, 'rCIP') || strcmp(CVI, 'NI')
    epsilon = 12;
    delta = 10^(-epsilon/dim);
    valind_batch.delta_term = delta.*eye(dim, dim);
end
fprintf('Done.\n');
fprintf('************************************************************************************\n')
fprintf('************************************************************************************\n')
  
%% Experimental example
fprintf('Clustering...\n');
rng(1, 'twister');  
options = statset('MaxIter', 1000);
k = 2:5;
len_k = length(k);
idx = zeros(nSamples, len_k);
CVI_value = zeros(1, len_k);
SFAM = cell(1, len_k);
maxEpochs = 1; % maximum number of epochs for SFAM (online = 1)
xcc = [data 1 - data];  % Complemenent-code 
CONN = cell(1, len_k);
for ix=1:len_k
    % kmeans
%     idx(:, ix) = kmeans(data, k(ix),'MaxIter',100, 'Replicates', 10, 'Start', 'plus'); 
    % EM
    gmfit = fitgmdist(data, k(ix), 'CovarianceType', 'full', 'SharedCovariance', false, 'Options', options);
    idx(:, ix) = cluster(gmfit, data); 
    % Compute CVI
    if ~strcmp(CVI, 'CONN') 
        valind_batch = valind_batch.param_batch(data, idx(:, ix));
    else      
        [cl, inds] = sort(idx(:, ix), 'ascend');
        xcc_sorted = xcc(inds,:);
        SFAM{ix} = SFuzzyARTMAPMod(settings); 
        SFAM{ix}.display = false; 
        SFAM{ix}.ART_A.find_bmu2 = true;
        for jx = 1:nSamples
            SFAM{ix} = SFAM{ix}.train(xcc_sorted(jx, :), maxEpochs, cl(jx)); 
        end
        [labels_A_1, labels_A_2] = SFAM{ix}.ART_A.classify(xcc);
        valind_batch = valind_batch.param_batch(labels_A_1, labels_A_2, idx(:, ix));
        CONN{ix} = valind_batch.CONN;
    end
    valind_batch = valind_batch.evaluate();  
    CVI_value(ix) = valind_batch.CriterionValue;
end
max_optimal = {'CH', 'SIL', 'I', 'GD43', 'GD53', 'rH', 'PS', 'CONN'};
min_optimal = {'WB', 'XB', 'DB', 'NI','rCIP'};
if ismember(CVI, max_optimal)
    [~, best_ind] = max(CVI_value);
elseif ismember(CVI, min_optimal)
    [~, best_ind] = min(CVI_value);    
end 
fprintf('Done.\n');

%% Figures
% Figure parameters
FONTSIZE = 16;
FONTWEIGHT = 'bold';
LINEWIDTH = 2;

fig = figure('visible', 'on');
hold on
% CVI
subplot(1,2,1)
stem(k, CVI_value, 'LineWidth', LINEWIDTH)
xlim([k(1)-1 k(end)+1])
title(CVI_name, 'Interpreter', 'none')
xlabel('Number of Clusters')
ylabel('CVI Value')
set(gca,'FontSize',FONTSIZE,'FontWeight',FONTWEIGHT)
box on

% Data
subplot(1,2,2)
title('Best Solution')
hold on   
box on 
if dim==2 
    clrs = rand(k(best_ind), 3);
    if strcmp(CVI, 'CONN')    
        % Prototypes & data for CONN_Index
        draw_FA_categories(SFAM{best_ind}.ART_A, data, idx(:, best_ind), 2, SFAM{best_ind}.map, false, clrs)
        hold on
        CONNvis(SFAM{best_ind}.ART_A, SFAM{best_ind}.map, CONN{best_ind}, clrs);
        set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]);  
    else            
        for ix=1:k(best_ind)
            plot(data(idx(:, best_ind)==ix,1),data(idx(:, best_ind)==ix,2), 'Marker', 'o', 'MarkerSize', 10, 'Color', clrs(ix,:), 'LineStyle', 'None');
        end
    end
end
% Axis properties
box on
set(gcf, 'Color','w'); 
set(gca, 'Color', 'none');
ax = gca;
ax.FontSize = FONTSIZE;
ax.FontWeight = FONTWEIGHT;

% Fullscreen
set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1])  