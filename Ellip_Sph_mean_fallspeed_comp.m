%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Article: "How Snow Aggregate Shapes and 
% Orientations Affects Fall Speed and Self-
%Collection Rates"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script: Ellip_Sph_mean_fallspeed_comp.m
%
% Description: This script calculates the bulk 
% fall speed quantities (number-, mass-, and reflectivity-weighted
% fall speeds) used in Figure 7 and 8 in the main article.
% See Table 1 in the main article for explicit forms
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% test_switch
% TEST #:
% 1: Disperse mass and projected area (Figure 4)
% 2: Same mass distribution (Figure 6)
% 3: Same projected area distribution (Figure 5)
test_switch = 1;

H_test = [2 0 3 1];

% Colors
b1 = [0.00 0.45 0.74];
r1 = [0.64 0.08 0.18];

g = 9.81;
Rd = 287.5;
P = 1e5;

% mks
TC = -5;
TK = 273.15+TC;
rhoa = P./(Rd.*TK);
eta_a = (1.496e-6 .* TK.^(3./2) ./(TK+120));

% Convert to cgs
g = g.*100;
rhoa = rhoa./1000;
eta_a = 10.*eta_a;

% AR_circ parameterization (i.e. Heymsfield 2002a,b)
n = 1.5;
alph = -0.6;

if alph == -0.8
    k = 0.0540;
elseif alph == -0.6
    k=0.078;
else
    k = 0.015;
    alph = -1;
end

% 'k' factor for area ratio
%kc = 1;
kc = [1 1 1 0.5];

% Area parameterization (i.e. Mitchell 1996) cgs
zet = 1.88;
sig = 0.2285;
% Size distribution parameters
nu = 1;

Lambda = logspace(log10(2),2,50);

Dn = 1./Lambda;
an = Dn./2;

% Aggregate projected area exponents
bet_ba = 1.0;
bet_ca = 0.;

if test_switch == 3
    bet_ba = 0;
    bet_ca = 0;
end

% MASC ellipsoid parameters
a_ba = 6.9793;
b_ba = 4.3502;
b_cb = 5.3437;

% mean shape
phimean = beta_moms(1,1,a_ba,b_ba,b_cb); % 1st moment
phivar = beta_moms(2,2,a_ba,b_ba,b_cb); % 2nd moment

bet_ar = zet-2;
alph_ar = 2.^zet .* sig./pi;

% SPHERES

% params for m-d relationship
bet_ms = 3+n.*bet_ar+alph;
%NOTE: alph_ms should probably have 10^p somewhere in here.
alph_ms = (2.^(2+n.*(zet)+alph)./3) .* pi.^(1-n) .* k .* sig.^n;

bm_ellip = NaN(length(H_test),length(an));
am_ellip = NaN(length(H_test),length(an));
bm_sph = NaN(length(H_test),length(an));
am_sph = NaN(length(H_test),length(an));

alph_v_sph = NaN(length(H_test),length(an));
alph_v_ellip = NaN(length(H_test),length(an));
bet_v_sph = NaN(length(H_test),length(an));
bet_v_ellip = NaN(length(H_test),length(an));

zet_ba = NaN(length(H_test),length(an));
zet_ca = NaN(length(H_test),length(an));

vtn_sph = NaN(length(H_test),length(an));
vtn_ellip = NaN(length(H_test),length(an));
vtm_sph = NaN(length(H_test),length(an));
vtm_ellip = NaN(length(H_test),length(an));
vtz_sph = NaN(length(H_test),length(an));
vtz_ellip = NaN(length(H_test),length(an));

X_n_sph = NaN(length(H_test),length(an));
X_n_ellip = NaN(length(H_test),length(an));
X_bar_sph = NaN(length(H_test),length(an));
X_bar_ellip = NaN(length(H_test),length(an));

alph_x = NaN(length(H_test),length(an));
alph_xs = NaN(length(H_test),length(an));

%if test_switch == 1 || test_switch==3
%alph_x = NaN(length(H_test),length(an));
%alph_xs = NaN(length(H_test),length(an));
%else
 %   alph_x = NaN(1,length(H_test));
  %  alph_xs = NaN(1,length(H_test));
%end

bet_x = NaN(1,length(H_test));
bet_xs = NaN(1,length(H_test));
chi_ba = NaN(1,length(H_test));
chi_ca = NaN(1,length(H_test));

if test_switch == 1 || test_switch == 3
    
mn_s = alph_ms.*an.^bet_ms;
    
q_bar = gamma(nu+bet_ms)./gamma(nu);

Z_bar = gamma(nu+2.*bet_ms)./gamma(nu);

fgam = @(b) gammaln(nu+6+2.*b)-2.*gammaln(nu+3+b)-log((Z_bar./q_bar.^2).*phimean.^2/phivar);

bet_rho_ellip = fzero(fgam,[-2 0]);
bet_m = bet_rho_ellip+3;
                           
elseif test_switch == 2
    
    alph_m = alph_ms.*ones(1,length(an));
    bet_m = bet_ms;
    
    mn_s = alph_m.*an.^bet_m;
    
end

for i = 1:length(H_test)

                
alph_rho_ellip =  mn_s.*q_bar./((4./3).*pi.*an.^(bet_m)).*...
                               (gamma(nu)./gamma(nu+bet_m)).*...
                               (1./phimean);
       
if test_switch == 1 || test_switch == 3                           
alph_m = (4./3).*pi.*alph_rho_ellip;        
end
    
    
alph_x(i,:) = (8./pi) .*g.*(rhoa./eta_a.^2).*alph_m.*...
                  alph_ar.^(-kc(i));
              
bet_x(i) = bet_m-kc(i).*bet_ar;             

bet_xs(i) = bet_ms-kc(i).*bet_ar;

alph_xs(i) = (8./pi) .*g.*(rhoa./eta_a.^2).*alph_ms.*...
                  alph_ar.^(-kc(i));
              
% Aggregate Best number ellipsoid coefficients
if test_switch == 1 || test_switch == 3
chi_ba(i) = kc(i).*bet_ba+1;
chi_ca(i) = kc(i).*bet_ca+1;
elseif test_switch == 2
chi_ba(i) = kc(i).*bet_ba;
chi_ca(i) = kc(i).*bet_ca;
end              

X_n_sph(i,:) = alph_xs(i).*an.^bet_xs(i);
X_n_ellip(i,:) = alph_x(i).*an.^bet_x(i);

% Mean Best number if using ellipsoid shell and circle area ratio param
X_bar_sph(i,:) = X_n_sph(i,:).*gamma(nu+bet_xs(i))./gamma(nu);
X_bar_ellip(i,:) = X_n_ellip(i,:).*gamma(nu+bet_x(i))./gamma(nu).*...
   beta_moms(chi_ba(i),chi_ca(i),a_ba,b_ba,b_cb);

% Solve for Reynold's number power-law coefficients

[bm_ellip(i,:),am_ellip(i,:)] = fallspeed_params(H_test(i),X_bar_ellip(i,:));
[bm_sph(i,:),am_sph(i,:)] = fallspeed_params(H_test(i),X_bar_sph(i,:));

%alph_vt for area ratio param.
alph_v_sph(i,:) = 0.5.*(eta_a./rhoa).*am_sph(i,:).*alph_xs(i).^bm_sph(i,:);
alph_v_ellip(i,:) = 0.5.*(eta_a./rhoa).*am_ellip(i,:).*alph_x(i).^bm_ellip(i,:);

bet_v_sph(i,:) = bm_sph(i,:).*bet_xs(i)-1;
bet_v_ellip(i,:) = bm_ellip(i,:).*bet_x(i)-1;

zet_ba(i,:) = chi_ba(i).*bm_ellip(i,:)-bet_ba./2;
zet_ca(i,:) = chi_ca(i).*bm_ellip(i,:)-bet_ca./2;

vt0n_ellip = alph_v_ellip(i,:).*an.^bet_v_ellip(i,:);
vt0n_sph = alph_v_sph(i,:).*an.^bet_v_sph(i,:);

vtn_sph(i,:) =0.01.* alph_v_sph(i,:) .*an.^(bet_v_sph(i,:)) .*(gamma(nu+bet_v_sph(i,:))./gamma(nu));

vtn_ellip(i,:) = 0.01.*alph_v_ellip(i,:) .*an.^(bet_v_ellip(i,:)) .*(gamma(nu+bet_v_ellip(i,:))./gamma(nu)) .*...
    beta_moms(zet_ba(i,:),zet_ca(i,:),a_ba,b_ba,b_cb);

vtm_sph(i,:) = 0.01.*alph_v_sph(i,:) .* an.^(bet_v_sph(i,:)) .* (gamma(nu+bet_v_sph(i,:)+bet_ms)./gamma(nu+bet_ms));

vtm_ellip(i,:) = 0.01.*alph_v_ellip(i,:).*an.^(bet_v_ellip(i,:)) .*(gamma(nu+bet_v_ellip(i,:)+bet_m)./gamma(nu+bet_m)) .*...
    beta_moms(zet_ba(i,:)+1,zet_ca(i,:)+1,a_ba,b_ba,b_cb)./phimean;

vtz_sph(i,:) = 0.01.*(alph_v_sph(i,:) .* an.^(bet_v_sph(i,:))) .* (gamma(nu+bet_v_sph(i,:)+2.*bet_ms)./gamma(nu+2.*bet_ms));

vtz_ellip(i,:) = 0.01.*(alph_v_ellip(i,:).*an.^(bet_v_ellip(i,:))) .*(gamma(nu+bet_v_ellip(i,:)+2.*bet_m)./gamma(nu+2.*bet_m)) .*...
    beta_moms(zet_ba(i,:)+2,zet_ca(i,:)+2,a_ba,b_ba,b_cb)./phivar;




end

% PLOT STUFF

figure;
subplot(3,1,1)
plot(Lambda,vtn_sph(1,:),'Color',r1,'linewidth',2.0);
hold on;
plot(Lambda,vtn_ellip(1,:),'Color',b1,'linewidth',2.0);
plot(Lambda,vtn_sph(2,:),'--','Color',r1,'linewidth',2.0);
plot(Lambda,vtn_ellip(2,:),'--','Color',b1,'linewidth',2.0);
plot(Lambda,vtn_sph(3,:),'-.','Color',r1,'linewidth',2.0);
plot(Lambda,vtn_ellip(3,:),'-.','Color',b1,'linewidth',2.0);
plot(Lambda,vtn_sph(4,:),':','Color',r1,'linewidth',2.0);
plot(Lambda,vtn_ellip(4,:),':','Color',b1,'linewidth',2.0);

set(gca,'xscale','log','xdir','reverse')
set(gca,'xscale','log')
xlim([2 100])
ylim([0 3.5])


subplot(3,1,2)
plot(Lambda,vtm_sph(1,:),'Color',r1,'linewidth',2.0);
hold on;
plot(Lambda,vtm_ellip(1,:),'Color',b1,'linewidth',2.0);
plot(Lambda,vtm_sph(2,:),'--','Color',r1,'linewidth',2.0);
plot(Lambda,vtm_ellip(2,:),'--','Color',b1,'linewidth',2.0);
plot(Lambda,vtm_sph(3,:),'-.','Color',r1,'linewidth',2.0);
plot(Lambda,vtm_ellip(3,:),'-.','Color',b1,'linewidth',2.0);
plot(Lambda,vtm_sph(4,:),':','Color',r1,'linewidth',2.0);
plot(Lambda,vtm_ellip(4,:),':','Color',b1,'linewidth',2.0);

set(gca,'xscale','log','xdir','reverse')
set(gca,'xscale','log')
xlim([2 100])
ylim([0 3.5])


subplot(3,1,3)
plot(Lambda,vtz_sph(1,:),'Color',r1,'linewidth',2.0);
hold on;
plot(Lambda,vtz_ellip(1,:),'Color',b1,'linewidth',2.0);
plot(Lambda,vtz_sph(2,:),'--','Color',r1,'linewidth',2.0);
plot(Lambda,vtz_ellip(2,:),'--','Color',b1,'linewidth',2.0);
plot(Lambda,vtz_sph(3,:),'-.','Color',r1,'linewidth',2.0);
plot(Lambda,vtz_ellip(3,:),'-.','Color',b1,'linewidth',2.0);
plot(Lambda,vtz_sph(4,:),':','Color',r1,'linewidth',2.0);
plot(Lambda,vtz_ellip(4,:),':','Color',b1,'linewidth',2.0);

set(gca,'xscale','log','xdir','reverse')
set(gca,'xscale','log')
xlim([2 100])

ylim([0 3.5])

figure;
plot(Lambda,vtm_sph(1,:)./vtn_sph(1,:),'Color',r1);
hold on;
plot(Lambda,vtm_ellip(1,:)./vtn_ellip(1,:),'Color',b1);
plot(Lambda,vtm_sph(2,:)./vtn_sph(2,:),'--','Color',r1);
plot(Lambda,vtm_ellip(2,:)./vtn_ellip(2,:),'--','Color',b1);
plot(Lambda,vtm_sph(3,:)./vtn_sph(3,:),'-.','Color',r1);
plot(Lambda,vtm_ellip(3,:)./vtn_ellip(3,:),'-.','Color',b1);
plot(Lambda,vtm_sph(4,:)./vtn_sph(4,:),':','Color',r1);
plot(Lambda,vtm_ellip(4,:)./vtn_ellip(4,:),':','Color',b1);

set(gca,'xscale','log','xdir','reverse')
set(gca,'xscale','log')
xlim([2 100])

figure;
plot(Lambda,vtm_ellip(1,:)./vtn_ellip(1,:).*(vtm_sph(1,:)./vtn_sph(1,:)).^(-1),'Color','k');
hold on;
plot(Lambda,vtm_ellip(2,:)./vtn_ellip(2,:).*(vtm_sph(2,:)./vtn_sph(2,:)).^(-1),'--','Color','k');
plot(Lambda,vtm_ellip(3,:)./vtn_ellip(3,:).*(vtm_sph(3,:)./vtn_sph(3,:)).^(-1),'-.','Color','k');
plot(Lambda,vtm_ellip(4,:)./vtn_ellip(4,:).*(vtm_sph(4,:)./vtn_sph(4,:)).^(-1),':','Color','k');

set(gca,'xscale','log','xdir','reverse')
set(gca,'xscale','log')
xlim([2 100])



function [bm,am] = fallspeed_params(Heyms_switch,Xbar)
% Aggregate fallspeed parameters 

%Heyms_switch = 0;

if Heyms_switch == 1
% Heymsfield and Westbrook 2010
ao = 0;
bo = 0;
Co = 0.35;
delta_0 = 8.0;
dbm = 0;
Psi_tX = 1.0;
elseif Heyms_switch == 2
% Mitchell1996
ao = 0;
bo = 0;
Co = 0.6;
delta_0 = 5.83;
dbm = 0;
Psi_tX = 1.0;
elseif Heyms_switch == 0
% Mitchell and Heymsfield 2005
ao = 1.7e-3;
bo = 0.8;
Co = 0.6;
delta_0 = 5.83;
dbm = 0;
Psi_tX = 1.0;
elseif Heyms_switch == 3
%KC2005
ao = 0;
bo = 0;
Co = 0.6;
delta_0 = 5.83;
X0 = 2.8e6;
k0 = 2;
Ct = 1.3;
z = Xbar./X0;

Psi_X = (1+z.^k0)./(1+Ct.*z.^k0);

dbm = - (k0.*(Ct-1).*z.^k0)./...
    (2.*(1+z.^k0).*(1+Ct.*z.^k0));

Psi_tX = sqrt(Psi_X)./(Xbar.^dbm);

end

C2 = (delta_0.^2) ./4;
C1 = 1./(C2.*sqrt(Co));

bml = (C1.*sqrt(Xbar)./...
(2.*(sqrt(1+C1.*sqrt(Xbar)) - 1).*sqrt(1+C1.*sqrt(Xbar)))) -...
((ao.*bo.*Xbar.^bo)./(C2.*(sqrt(1+C1.*sqrt(Xbar))-1).^2));

bm = bml+dbm;

am = ...
Psi_tX.*(C2.*((sqrt(1+C1.*sqrt(Xbar))-1).^2) - ao.*Xbar.^bo)./(Xbar.^(bml));

end


