%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Article: "How Snow Aggregate Shapes and 
% Orientations Affects Fall Speed and Self-
%Collection Rates"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fallspeed_agg_MC_test_ellipproj.m
% Description: Calculates Fall speed distribution
% used in Figures 4-6 of main article.

% Changed to use rho_ellip instead of rho_sph

a_switch = 0;

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

% Size distribution parameters
nu = 1;

% MASC ellipsoid parameters
a_ba = 6.9793;
b_ba = 4.3502;
b_cb = 5.3437;

% mean shape
phiba_phica_bar = beta_moms(1,1,a_ba,b_ba,b_cb);

% AR_circ parameterization (i.e. Heymsfield 2002a,b)
%k = 0.01;
%n = 1.5;
n = 1.52;
k = 0.18;
alph = -0.8;

% Area parameterization (i.e. Mitchell 1996) cgs
zet = 1.88;
sig = 0.2285;

% params for m-d relationship
bet_m = 3+n.*(zet-2)+alph;
alph_m = (2.^(2+n.*(zet)+alph)./3) .* pi.^(1-n) .* k .* sig.^n;

% If same mass
q_agg = 1e-3 ;
N_agg = 0.001;
q_bar = q_agg/N_agg;

%na = makedist('gamma','a',nu,'b',an);
if a_switch == 0
a_bar = 0.0238*4; % centimeters
%a_bar = 1;

an = a_bar./nu;
an_sph = an;
an_ellip = an;

elseif a_switch == 1
    an_sph = ((rhoa.*q_bar./alph_m).*gamma(nu)./...
    gamma(nu+bet_m)).^(1./bet_m);

    an_ellip = an_sph .* phiba_phica_bar.^(-1./bet_m);
     
end

na_sph = truncate(makedist('gamma','a',nu,'b',an_sph),0,100.*an_sph);
na_ellip = truncate(makedist('gamma','a',nu,'b',an_ellip),0,100.*an_ellip);


% params for x-d relationship
bet_x = bet_m+(2-zet)./4;
alph_x = 8.*alph_m.*g.*rhoa./(pi.*eta_a.^2) .*2.^(-zet./4) .*...
sig.^(-1./4) .*pi.^(1./4);

X_n_sph = alph_x.*an_sph.^bet_x;
X_n_ellip = alph_x.*an_ellip.^bet_x;

% Aggregate projected area exponents
bet_ba = 1;
bet_ca = 0;

% Aggregate Best number ellipsoid coefficients


% Mean Best number if using ellipsoid shell and circle area ratio param
X_bar_sph = X_n_sph.*gamma(nu+bet_x)./gamma(nu);
%X_bar_ellip = X_n_ellip.*gamma(nu+bet_x)./gamma(nu).*...
 %  beta_moms(1+bet_ba./4,1+bet_ca./4,a_ba,b_ba,b_cb);

 X_bar_ellip = X_n_ellip.*gamma(nu+bet_x)./gamma(nu).*...
   beta_moms(1+bet_ba.*(0.25-n),1+bet_ca.*(0.25-n),a_ba,b_ba,b_cb);
 

% Solve for Reynold's number power-law coefficients
[bm_ellip,am_ellip] = fallspeed_params(X_bar_ellip);
[bm_sph,am_sph] = fallspeed_params(X_bar_sph);

%alph_vt for area ratio param.
alph_v_sph = 0.5.*(eta_a./rhoa).*am_sph.*alph_x.^bm_sph;
alph_v_ellip = 0.5.*(eta_a./rhoa).*am_ellip.*alph_x.^bm_ellip;

bet_v_sph = bm_sph.*bet_x-1;
bet_v_ellip = bm_ellip.*bet_x-1;

[nphi_biv,phib_bins,phic_bins] = nphi_biv_agg(a_ba,b_ba,b_cb);

nphi_biv(isnan(nphi_biv)) = 0;

phiba_space = 0.05:0.05:1.0;
phica_space = phiba_space;

N_biv_total = 20000;
phiba_samp = NaN(1,N_biv_total);
phica_samp = phiba_samp;

for i = 1 : N_biv_total
    
 [phiba_samp(i),phica_samp(i)] = pinky(phib_bins,phic_bins,nphi_biv');   
    
end

a_sph_samp = random(na_sph,[1 N_biv_total]);
a_ellip_samp = random(na_ellip,[1 N_biv_total]);

AR_circ = 2.^(zet) .* pi.^(-1) .* sig .* a_sph_samp.^(zet-2);
AR_ellip = AR_circ./(phiba_samp.^(bet_ba).*phica_samp.^(bet_ca));

rhoi = k.*AR_circ.^(n+alph);

m_sph = alph_m.*a_sph_samp.^bet_m;
%m_ellip = alph_m.*a_ellip_samp.^bet_m .*phiba_samp.*phica_samp;
m_ellip = alph_m.*a_ellip_samp.^bet_m .*...
    phiba_samp.^(1-n.*bet_ba).*...
    phica_samp.^(1-n.*bet_ca);


vt_sph = alph_v_sph.*a_sph_samp.^bet_v_sph;

%zet_ba = bm_ellip+0.25.*bet_ba.*(bm_ellip-2);
%zet_ca = bm_ellip+0.25.*bet_ca.*(bm_ellip-2);

zet_ba = bm_ellip+(0.25-n).*bet_ba.*(bm_ellip-2);
zet_ca = bm_ellip+(0.25-n).*bet_ca.*(bm_ellip-2);


vt_ellip = alph_v_ellip.*a_ellip_samp.^(bet_v_ellip) .*...
    phiba_samp.^(zet_ba) .*...
    phica_samp.^(zet_ca);

vt0n_ellip = alph_v_ellip.*an_ellip.^bet_v_ellip;
vt0n_sph = alph_v_sph.*an_sph.^bet_v_ellip;

vt_min = min(vt_ellip);
vt_max = max(vt_ellip);

vt_norm = logspace(log10(vt_min),log10(vt_max),200);

% H-function representation
n_vt = (gamma(a_ba+b_ba+b_cb)./gamma(a_ba)) .* (1./gamma(nu)).*...
    Fox_H([nu-bet_v_ellip a_ba+b_ba-zet_ca a_ba-(zet_ba+zet_ca)],[bet_v_ellip zet_ca zet_ba+zet_ca],[],[],[a_ba+b_ba+b_cb-zet_ca a_ba+b_ba-(zet_ba+zet_ca)],[zet_ca zet_ba+zet_ca],[],[],vt_norm./vt0n_ellip);

%vtn_sph = alph_v_sph .*an_sph.^(bet_v_sph) .*(gamma(nu+bet_v_sph)./gamma(nu));

%vtn_ellip = alph_v_ellip .*an_ellip.^(bet_v_ellip) .*(gamma(nu+bet_v_ellip)./gamma(nu)) .*...
%    beta_moms(zet_ba,zet_ca,a_ba,b_ba,b_cb);

%vtm_sph = alph_v_sph .* an_sph.^(bet_v_sph) .* (gamma(nu+bet_v_sph+bet_m)./gamma(nu+bet_m));

%vtm_ellip = alph_v_ellip .*an_ellip.^(bet_v_ellip) .*(gamma(nu+bet_v_ellip+bet_m)./gamma(nu+bet_m)) .*...
%    beta_moms(zet_ba+1,zet_ca+1,a_ba,b_ba,b_cb)./phiba_phica_bar;

%msph_min = min(m_samp./mn);
%msph_max = max(m_samp./mn);

% Now use Fox H-function

%m_norm = logspace(log10(msph_min),log10(msph_max),200);

%m_H = (gamma(aba+bba+bcb)./gamma(aba)) .*(1./gamma(nu)) .*...
%    Fox_H([nu-bet_m aba-2 aba+bba-1],[bet_m 2 1],[],[],[aba+bba-2 aba+bba+bcb-1],[2 1],[],[],m_norm);

%figure;
%histogram(log10(m_samp./mn),'normalization','pdf');

%hold on;
%plot(log10(m_norm),m_norm.*m_H.*log(10));


figure;
histogram(log10(vt_sph./100),'normalization','pdf');
hold on;
histogram(log10(vt_ellip./100),'normalization','pdf');
hold on;
plot(log10(vt_norm./100),(vt_norm./vt0n_ellip).*n_vt.*log(10));
set(gca,'xtick',[-1 0 log10(2)],'xticklabel',{'0.1' '1.0' '2.0'})
xlim([-1 log10(2)])

figure;
histogram(log10(1000.*m_sph.*vt_sph./100),'normalization','pdf');
hold on;
histogram(log10(1000.*m_ellip.*vt_ellip./100),'normalization','pdf');
%set(gca,'xtick',[-3 -2 -1 0 1],'xticklabel',{'0.001' '0.01' '0.1' '1.0' '10.0'})
%xlim([-3 1])


%disp(vtn_sph./100)
%disp(vtn_ellip./100)
%disp(vtm_sph./100)
%disp(vtm_ellip./100)

%disp('--------------')
%disp(vtm_sph./vtn_sph)
%disp(vtm_ellip./vtn_ellip)
%disp((vtm_ellip./vtn_ellip)./(vtm_sph./vtn_sph))

%-----------------------------------------------------%

% This code does random sampling from the generalized gamma distribution
% P(R) = Beta*R^(alpha-1)*exp(-(R/a)^Beta)/(a^alpha*gamma(alpha/beta))

%function R = gengamrand(n,Alpha, Beta, a)

%P = rand(n);
%for i = 1:length(Alpha)
%    for j = 1:length(a)
%        X(:,i,j)= gammaincinv(P,(Alpha(i)/Beta));
%        R(:,i,j) = a(j).* (X(:,i,j)).^(1/Beta);
%    end
%end
%end

%-----------------------------------------------------%

function [bm,am] = fallspeed_params(Xbar)
% Aggregate fallspeed parameters 
% Mitchell and Heymsfield 2005
ao = 1.7e-3;
bo = 0.8;
Co = 0.6;
delta_0 = 5.83;

C2 = (delta_0.^2) ./4;
C1 = 1./(C2.*sqrt(Co));

bm = ...
(C1.*sqrt(Xbar)./...
(2.*(sqrt(1+C1.*sqrt(Xbar)) - 1).*sqrt(1+C1.*sqrt(Xbar)))) -...
((ao.*bo.*Xbar.^bo)./(C2.*(sqrt(1+C1.*sqrt(Xbar))-1).^2));

am = ...
(C2.*((sqrt(1+C1.*sqrt(Xbar))-1).^2) - ao.*Xbar.^bo)./(Xbar.^bm);


end

