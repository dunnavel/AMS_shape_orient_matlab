%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Article: "How Snow Aggregate Shapes and 
% Orientations Affects Fall Speed and Self-
%Collection Rates"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proj_plaw.m
%
% Description: Estimates geometric factors for 
% average projected area (Psi) and projected 
% semi-major dimension (Phi) used in Figure
% 1 of main article.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100000;

% phi space
phi_samp = linspace(0.1,1,phimax);

% power-law estimates for "random" orientations
psi_ba = 0.9;
psi_ca = 0.5;
phi_ba = 0.3;
phi_ca = 0.03;

theta_rand = 2.*pi.*rand(1,N);

phi_rand = acosd(2*rand(1,N)-1);

phimax = 100;

l2 = @(theta,phi) cos(theta).^2 .* sin(phi).^2;
m2 = @(theta,phi) sin(theta).^2 .* sin(phi).^2;
n2 = @(theta,phi) cos(phi).^2;

l2_samp = l2(theta_rand,phi_rand);
m2_samp = m2(theta_rand,phi_rand);
n2_samp = n2(theta_rand,phi_rand);

Psi_proj = @(l2,m2,n2,phiba,phica) ...
    sqrt(l2.*phica.^2+m2.*phiba.^2+n2.*phiba.^2.*phica.^2);

phi_proj = @(l2,m2,n2,phiba,phica) ...
 l2.*(1+phica.^2)+m2.*(1+phiba.^2)+n2.*(phiba.^2+phica.^2);   

Phi_L_proj = @(l2,m2,n2,phiba,phica) ...
 sqrt((phi_proj(l2,m2,n2,phiba,phica) + ...
 sqrt(phi_proj(l2,m2,n2,phiba,phica).^2 - 4.* Psi_proj(l2,m2,n2,phiba,phica)))./2);

 Phi_test_proj = @(l2,m2,n2,phiba,phica) ...
     phi_proj(l2,m2,n2,phiba,phic);
  
Psi_test = NaN(phimax,phimax);
Phi_out = NaN(phimax,phimax);
Phi_test = NaN(phimax,phimax);

Phi_plaw = NaN(phimax,phimax);
Psi_plaw = NaN(phimax,phimax);

for i = 1:phimax
    for j = 1:phimax
        
        if i>=j
            
            Psi_test(i,j) = ...
            mean(Psi_proj(l2_samp,m2_samp,n2_samp,phi_samp(i),phi_samp(j)),'omitnan');
               
           % Phi_out(i,j) = ...
           % mean(Phi_L_proj(l2_samp,m2_samp,n2_samp,phi_samp(i),phi_samp(j)),'omitnan');
                
            phi_temp = phi_proj(l2_samp,m2_samp,n2_samp,phi_samp(i),phi_samp(j));
        
            phi_L_temp = sqrt((phi_temp + sqrt(phi_temp.^2 - 4.*Psi_proj(l2_samp,m2_samp,n2_samp,phi_samp(i),phi_samp(j)).^2))./2);
            
            Phi_test(i,j) = mean(phi_L_temp,'omitnan');
            
            Phi_plaw(i,j) = phi_samp(i).^(phi_ba) .* phi_samp(j).^(phi_ca);
            Psi_plaw(i,j) = phi_samp(i).^(psi_ba).*phi_samp(j).^(psi_ca);
            
            
        end
        
        
    end
    
end

Psi_test =Psi_test';
Phi_test = real(Phi_test');

figure; 
subplot(1,2,1);
contourf(phi_samp,phi_samp,Phi_test);
caxis([0.5 1.0])
subplot(1,2,2);
contourf(phi_samp,phi_samp,Phi_plaw');
caxis([0.5 1.0])

figure; 
subplot(1,2,1);
contourf(phi_samp,phi_samp,Psi_test);
caxis([0.1 1.0])
subplot(1,2,2);
contourf(phi_samp,phi_samp,Psi_plaw');
caxis([0.1 1.0])

figure; 
subplot(1,2,1);
contourf(phi_samp,phi_samp,Phi_test-Phi_plaw');
subplot(1,2,2);
contourf(phi_samp,phi_samp,Psi_test-Psi_plaw');


