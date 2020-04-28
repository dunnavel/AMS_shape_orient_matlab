%betm_gam_plaw.m

% MASC ellipsoid parameters
a_ba = 6.9793;
b_ba = 4.3502;
b_cb = 5.3437;

nu = 1;

% Colors
b1 = [0.00 0.45 0.74];
r1 = [0.64 0.08 0.18];

% AR_circ parameterization (i.e. Heymsfield 2002a,b)
n = 1.5;
alph = -1;

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

bet_ar = zet-2;
alph_ar = 2.^zet .* sig./pi;

kc = 1.0;

bet_ba = 1.0;
bet_ca = 0;

% mean shape
phimean = beta_moms(1,1,aba,bba,bcb); % 1st moment
phivar = beta_moms(2,2,aba,bba,bcb); % 2nd moment

%bet_ms = linspace(1.4,2.6,50);
%bet_ms = [1.8 2.0 2.2 2.4 2.6 2.8 3.0];
bet_ms = [1.8 2.0 2.2 3.0];
%bet_ms = 3.0;
bet_m = nan(1,length(bet_ms));

for bi = 1 : length(bet_ms)
    
q_bar =  gamma(nu+bet_ms(bi))./gamma(nu);

Z_bar = gamma(nu+2.*bet_ms(bi))./gamma(nu);

fgam = @(b) gammaln(nu+6+2.*b)-2.*gammaln(nu+3+b)-log((Z_bar./q_bar.^2).*phimean.^2/phivar);


bet_rho_ellip = fzero(fgam,[-2 0]);
bet_m(bi) = bet_rho_ellip+3;

end

%bm_new = linspace(0.4,0.6,3);
%bet_vs = linspace(0,1.0,100);

bet_vs = logspace(-3,0,100);

%bet_vs = 0.3;

bm_new = NaN(length(bet_vs),length(bet_ms));
bet_ve = bet_vs;
ss_e = bm_new;
ss_s = ss_e;
%phi_factor = bm_new;

for i = 1 :length(bet_vs)
    for j = 1 : length(bet_ms)
        %bet_vs(i,j) = (bet_ms(j)-bet_ar).*bm_new(i)-1;
        bm_new(i,j) = (bet_vs(i)+1)./(bet_ms(j)-bet_ar);
        %phi_factor(i,j) = beta_moms(2.*bm_new(i,j)-0.5,bm_new(i,j),aba,bba,bcb);

        bet_ve(i,j) = (bet_m(j)-bet_ar).*bm_new(i,j)-1;
        
        zet_ba = (kc.*bet_ba+1).*bm_new(i,j)-bet_ba./2;
        zet_ca= (kc.*bet_ca+1).*bm_new(i,j)-bet_ca./2;
        
        %phi_factor = beta_moms(2.*zet_ba+2,2.*zet_ca+2,a_ba,b_ba,b_cb)./...
        %                       (beta_moms(zet_ba+1,zet_ca+1,a_ba,b_ba,b_cb)).^2;
        
        phi_factor = beta_moms(3.*zet_ba+1,3.*zet_ca+1,a_ba,b_ba,b_cb).*...
                              beta_moms(zet_ba+1,zet_ca+1,a_ba,b_ba,b_cb)./...
                              (beta_moms(2.*zet_ba+1,2.*zet_ca+1,a_ba,b_ba,b_cb).^2);
                           
       % gam_factor_ellip = gamma(nu+2.*(bet_m(j)+bet_ve(i,j))).*gamma(nu+bet_m(j)).^2 ./...
        %    (gamma(nu+2.*bet_m(j)).*gamma(nu+bet_ve(i,j)+bet_m(j)).^2);
        
        % gam_factor_sph = gamma(nu+2.*(bet_ms(j)+bet_vs(i))).*gamma(nu+bet_ms(j)).^2 ./...
       %     (gamma(nu+2.*bet_ms(j)).*gamma(nu+bet_vs(i)+bet_ms(j)).^2);
                           
       
       gam_factor_ellip = gamma(nu+bet_m(j)+3.*bet_ve(i,j)).*...
                                        gamma(nu+bet_m(j)+bet_ve(i,j))./...
                                        (gamma(nu+bet_m(j)+2.*bet_ve(i,j)).^2);
                                    
        gam_factor_sph = gamma(nu+bet_ms(j)+3.*bet_vs(i)).*...
                                        gamma(nu+bet_ms(j)+bet_vs(i))./...
                                        (gamma(nu+bet_ms(j)+2.*bet_vs(i)).^2);
       
        ss_e(i,j) = sqrt(gam_factor_ellip.*phi_factor-1);
        ss_s(i,j)=  sqrt(gam_factor_sph-1);
        
        %ss_e(i,j) = phi_factor(i,j).*gamma(1+bet_m(j)+bet_ve(i,j))./gamma(1+bet_ve(i,j));
        %ss_s(i,j) = gamma(1+bet_ms(j)+bet_vs(i))./gamma(1+bet_vs(i));
    end
end

ss_norm = ss_e./ss_s;

figure;
subplot(2,1,1)
plot(bet_vs,ss_s,'Color',r1);
hold on;
plot(bet_vs,ss_e,'Color',b1);
subplot(2,1,2);
plot(bet_vs,ss_norm,'Color','k');

% subplot(1,2,1);
% contourf(bet_ms,bet_vs,log10(ss_s));
% caxis(log10([0.001 6]))
% subplot(1,2,2);
% contourf(bet_ms,bet_vs,log10(ss_e));
% caxis(log10([0.001 6]))


%figure;
%plot(bet_ms,bet_vs,'b');
%hold on;
%plot(bet_ms,bet_ve,'g');




% params for m-d relationship
bet_ms = 3+n.*bet_ar+alph;

q_bar =  gamma(nu+bet_ms)./gamma(nu);

Z_bar = gamma(nu+2.*bet_ms)./gamma(nu);

fgam = @(b) gammaln(nu+6+2.*b)-2.*gammaln(nu+3+b)-log((Z_bar./q_bar.^2).*phimean.^2/phivar);


bet_rho_ellip = fzero(fgam,[-2 0]);
bet_m = bet_rho_ellip+3;



aref = 7.5;
a = linspace(0.05,10,50);
%a = 0.03:0.03:10;

err_b = 0.2.*ones(1,length(a));

bet_vs = 0.21;
alph_vs = 0.69;
%alph_vs = 0.7981;
%alph_vs = 0.4921;

bet_x = bet_m-bet_ar;

bet_xs = bet_ms-bet_ar;

bm = (bet_vs+1)./(bet_ms-bet_ar);

bet_ve = bet_x.*bm-1;

phi_factor = beta_moms(2.*bm-0.5,bm,aba,bba,bcb);

phi_2_factor = beta_moms(2.*(2.*bm-0.5),2.*bm,aba,bba,bcb);

phi_std = sqrt(phi_2_factor-phi_factor.^2);

alph_ve = alph_vs.*aref.^(bet_vs-bet_ve).*(1./phi_factor);

vt_s =alph_vs.*a.^bet_vs;

vt_e = alph_ve.*a.^bet_ve .* phi_factor;

vt_e_t = vt_e+phi_std;
vt_e_b = vt_e-phi_std;

%figure;
%plot(bet_ms,bet_m);

% figure; 
% errorbar(a,vt_s,err_b,'Color',r1,'linewidth',3.0);
% hold on;
% fill([a fliplr(a)],[vt_e_b fliplr(vt_e_t)],b1,'FaceAlpha',0.3,'EdgeColor','none');
% plot(a,vt_e,'Color',b1,'linewidth',3.0);
% hold on;
% ylim([0 1.5])
% set(gca,'ticklabelinterpreter','latex','fontsize',26,'layer','top','linewidth',1.4);
% set(gca,'xminortick','on','yminortick','on','ticklength',[0.025 0.05])
% set(gca,'ytick',0:0.2:1.5)
% grid on
% ylabel('Fall Speed (m/s)','interpreter','latex','fontsize',26)
% xlabel('Maximum Dimension (mm)','interpreter','latex','fontsize',26)
