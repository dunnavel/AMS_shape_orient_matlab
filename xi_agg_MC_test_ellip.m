%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Article: "How Snow Aggregate Shapes and 
% Orientations Affects Fall Speed and Self-
%Collection Rates"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xi_agg_MC_test_ellip.m
% Test out distribution of random variate
% xi = z^zeta_a * phiba^zeta_ba * phica^zeta_ca
% for ellipsoid aggregates using Monte Carlo
% simulation and H-function representation

% Total samples for MC
N_biv_total = 10000;

% Tswitch = 0, No size distribution truncation
% Tswitch = 1, Size distribution truncation
Tswitch = 0;

% Colors
b1 = [0.00 0.45 0.74];
r1 = [0.64 0.08 0.18];

% MASC ellipsoid parameters
aba = 6.9793;
bba = 4.3502;
bcb = 5.3437;

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


% Area parameterization (i.e. Mitchell 1996) cgs
zet = 1.88;
sig = 0.2285;

bet_ba = 1.0;
bet_ca =0.;

% a-axis size dist params
Ni = 1;
nu = 1;

Lambda = 10.0; %cgs
Dn = 1./Lambda;
ani = Dn./2;

% Swap in whatever power-law parameters you want.
zeta_a = zet-2;
zeta_ba = -bet_ba;
zeta_ca = -bet_ca;


    an = [];
    An = [];
    
    bm =  [nu-zeta_a aba-(zeta_ba+zeta_ca) aba+bba-zeta_ca];
    Bm =   [zeta_a zeta_ba+zeta_ca zeta_ca];
    
    ap = [aba+bba-(zeta_ba+zeta_ca) aba+bba+bcb-zeta_ca] ;
    Ap = [zeta_ba+zeta_ca zeta_ca]; 
    
    bq = [];
    Bq = [];


    % Switch H-function parameters
if zeta_a<=0 && zeta_ba<=0 && zeta_ca<=0
    
    an_new = 1-bm;
    An_new = abs(Bm);
    
    bq_new = 1- ap;
    Bq_new = abs(Ap);

    an = an_new;
    An = An_new;
    
    bm = [];
    Bm = [];
    
    ap = [];
    Ap =  [];
    
    bq = bq_new;
    Bq = Bq_new;

end

% Make sure the prefactor is consistent with 
% the characeteristic power law.
xi_pre = ...
10.^(2.*zet-4).* 2.^zet .*sig./pi ;

% An for projected area.
xin = xi_pre*ani.^(zeta_a);

% make size distribution
if Tswitch==0
% Un-truncated size distribution
na = truncate(makedist('gamma','a',nu,'b',ani),ani./1000,1000.*ani);
elseif Tswitch==1
% Truncated size distribution
na = truncate(makedist('gamma','a',nu,'b',ani),0.174/2,0.476/2);
end
a_samp = random(na,[1 N_biv_total]);

% Get aspect ratio distribution

[nphi_biv,phib_bins,phic_bins] = nphi_biv_agg(aba,bba,bcb);

nphi_biv(isnan(nphi_biv)) = 0;

phiba_space = 0.05:0.05:1.0;
phica_space = phiba_space;

phiba_samp = NaN(1,N_biv_total);
phica_samp = phiba_samp;

for i = 1 : N_biv_total
    
    
 [phiba_samp(i),phica_samp(i)] = pinky(phib_bins,phic_bins,nphi_biv');   
    
    
end

xi_samp = xi_pre*a_samp.^(zeta_a) .*...
                   phiba_samp.^(zeta_ba) .* ...
                   phica_samp.^(zeta_ca);

xi_min = min(xi_samp./xin);
xi_max = max(xi_samp./xin);

% Now use Fox H-function

%xi_norm = linspace(xi_min,xi_max,200);

xi_norm = logspace(log10(xi_min),log10(xi_max),200);

xi_H = (gamma((aba+bba+bcb))./gamma(aba)) .*(1./gamma(nu)) .*...
            (1./xin).* Fox_H(bm,Bm,an,An,ap,Ap,bq,Bq,xi_norm);

figure;
%histogram(log10(xi_samp),'normalization','pdf');
histogram(xi_samp,'normalization','pdf','Facecolor',r1);

hold on;
%plot(log10(xin.*xi_norm),xin.*xi_norm.*xi_H.*log(10));

plot(xin.*xi_norm,xi_H,'Color',r1,'linewidth',4.0);

