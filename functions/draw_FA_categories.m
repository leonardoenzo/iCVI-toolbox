%% """ Visualization of fuzzy ART hyperboxes """
% 
% PROGRAM DESCRIPTION
% This program plots fuzzy ART categories (hyperboxes), data samples, and 
% category labels according to the partition provided.
%
% INPUTS
% ART: trained SFAM class object
% data: data set matrix (rows: samples, columns: features)
% labels: partition labels
% lw: line width of hyperboxes
% Cat_map: category to class many-to-one binary mapping matrix (rows: categories, columns: classes)
% text_flag: true/false (display categories' indices)
% clrs: colors to be assigned to each class (matrix: number of classes by 3)
%
% REFERENCES
% [1] L. E. Brito da Silva, N. M. Melton, and D. C. Wunsch II, "Incremental
% Cluster Validity Indices for Hard Partitions: Extensions  and  Comparative  
% Study," ArXiv  e-prints, Feb 2019, arXiv:1902.06711v1 [cs.LG].
% 
% [2] G. Carpenter, S. Grossberg, and D. Rosen, "Fuzzy ART: Fast 
% stable learning and categorization of analog patterns by an adaptive 
% resonance system," Neural Networks, vol. 4, no. 6, pp. 759–771, 1991.
%
% Code written by Leonardo Enzo Brito da Silva
% Under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot categories and data
function draw_FA_categories(ART, data, labels, lw, Cat_map, text_flag, clrs)

    [~, dim] = size(data);
    ux = unique(labels); 
    nClasses = length(ux);
    
    % Comment this for loop if you have Statistics Toolbox and replace it
    % with gscatter built-in function
    for ix=1:nClasses
        plot(data(labels==ux(ix), 1), data(labels==ux(ix), 2),'Color', clrs(ix, :),...
            'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', clrs(ix, :),...
            'LineStyle', 'none');
    end  
    % Uncomment if you have Statistics Toolbox
    % gscatter(data(:, 1), data(:, 2), labels, clrs, '.', 10, 'off') 

    for j=1:ART.nCategories
        cat_color = clrs(Cat_map(j, :)==1, :);       

        if text_flag
            txt1 = num2str(j);    
            pos = ART.W(j,1:2) + (1-ART.W(j,3:4) - ART.W(j,1:2))./2;            
            text('Position',[pos(1), pos(2), 0],'String', txt1, 'BackgroundColor',cat_color,'FontSize',50,'FontWeight','bold','Color',[1 1 1])
        end
        
        if dim==2
            x = ART.W(j, 1);
            y = ART.W(j, 2);
            w = 1 - ART.W(j, 3) - ART.W(j, 1);
            h = 1 - ART.W(j, 4) - ART.W(j, 2);        
            if and((w>0), (h>0))
                pos = [x y w h]; 
                r = rectangle('Position', pos);
                r.FaceColor = 'none';
                r.EdgeColor = cat_color;
                r.LineWidth = lw;
                r.LineStyle = '-';
                r.Curvature = [0 0]; 
            else
                X = [ART.W(j, 1) 1 - ART.W(j, 3)];
                Y = [ART.W(j, 2) 1 - ART.W(j, 4)];
                l = line(X, Y);
                l.Color = cat_color;
                l.LineStyle = '-';
                l.LineWidth = lw;
                l.Marker = 'none';
            end
        end
    end

    axis square
    box on

    set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(gcf,'color','w');

end