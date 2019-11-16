% PROGRAM DESCRIPTION
% This is a MATLAB implementation of a similarity matrix graph visualizer. 
% If the input is the CONN matrix, then it becomes a CONNvis variant.
%
% REFERENCES:
% [1] K. Tasdemir and E. Merényi, "Exploiting data topology in visualization
% and  clustering  of  self-organizing  maps," IEEE  Trans.  Neural  Netw.,
% vol. 20, no. 4, pp. 549–562,  Apr. 2009.
%
% Code written by Leonardo Enzo Brito da Silva
% Under the supervision of Dr. Donald C. Wunsch II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function []  = CONNvis(ART, mapping, Mat, clrs)

    C = (ART.W(:,1:ART.dim) + (1 - ART.W(:,ART.dim+1:2*ART.dim)))./2;
    
    palette = colormap(flipud(copper)); 
    num_clrs_p = size(palette, 1);
    indices = find(Mat>0);
    
    M_color = Mat;  
    vals_color = Mat(indices);
    vals_color = mapminmax(vals_color', 1 , num_clrs_p);
    u_values = unique(vals_color);
    if length(u_values)==1
        vals_color = num_clrs_p.*ones(size(vals_color, 1), size(vals_color, 2));
    else
        vals_color = floor(vals_color);
    end 
    M_color(indices) = vals_color; 
    
    M_thickness = Mat;      
    vals_thickness = Mat(indices);
    vals_thickness = mapminmax(vals_thickness', 1 , 8);
    u_values = unique(vals_thickness);
    if length(u_values)==1
        vals_thickness = 8.*ones(size(vals_thickness, 1), size(vals_thickness, 2));
    end     
    M_thickness(indices) = vals_thickness;         

    for ix=1:ART.nCategories-1
        for jx=ix+1:ART.nCategories
            if Mat(ix, jx) > 0               
                plot(C([ix jx],1), C([ix jx],2),'Color',palette(M_color(ix, jx),:), 'Linewidth', M_thickness(ix, jx))     
            end
        end
    end
    
    h1 = colorbar('location','EastOutside','Ticks',[0 0.5 1], 'TickLabels',...
        {'Weak','Medium','Strong'},'FontWeight','bold','FontSize',20);
    
    nColors = size(mapping, 2);     
    for jx=1:ART.nCategories
        cat_color = clrs(mapping(jx, :)==1, :); 
        plot(C(jx, 1), C(jx, 2), 'Marker', 'p','MarkerFaceColor', cat_color,'MarkerSize',8, 'Color', cat_color)
    end
    
    axis square
    
end