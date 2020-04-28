%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Article: "How Snow Aggregate Shapes and 
% Orientations Affects Fall Speed and Self-
%Collection Rates"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%
% nphi_biv_agg.m
%
% Description:
% Calculates the Bivariate Beta distribution for 
% aggregate aspect ratios given in Equation 2
% in main paper.
%%%%%%%%%%%%%%%%%



function [nphi,phib_bins,phic_bins] = nphi_biv_agg(a_ba,b_ba,b_cb)

phib_bins = linspace(0.05,1.0,200);
phic_bins = linspace(0.05,1.0,200);


nphi = NaN(length(phib_bins),length(phic_bins));


nphi_func = @(phib,phic,a_ba,b_ba,b_cb) ...
    (1./beta(a_ba,b_ba)).*(1./beta(a_ba+b_ba,b_cb)) .* ...
    phic.^(a_ba+b_ba-1) .* (phib-phic).^(b_cb-1) .* ...
    phib.^(-b_ba-b_cb) .* (1-phib).^(b_ba-1);

for i = 1 : length(phib_bins)
    for j = 1 : length(phic_bins)
        
        
        
       if phic_bins(j) <= phib_bins(i)
           
            nphi(i,j) = nphi_func(phib_bins(i),phic_bins(j),...
                                  a_ba,b_ba,b_cb);
           
       end
        
        
    end
    
end










end