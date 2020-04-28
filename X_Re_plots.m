% X_Re_plots.m

fsize = 14;

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

nua =eta_a./rhoa;

kc = 1.0;


% params for m-d relationship
bet_ms = 1.8;
%NOTE: alph_ms should probably have 10^p somewhere in here.
alph_ms = 0.00145 ;

% Area parameterization (i.e. Mitchell 1996) cgs
zet = 1.88;
sig = 0.2285;

bet_ar = zet-2;
alph_ar = sig;

D = logspace(-2,1,200);

X = logspace(1,8,200);


[bm_0,am_0] = fallspeed_params(0,X);

[bm_1,am_1] = fallspeed_params(1,X);

[bm_2,am_2] = fallspeed_params(2,X);

[bm_3,am_3] = fallspeed_params(3,X);

bet_xs = bet_ms-kc.*bet_ar;

alph_xs = (2) .*g.*(rhoa./eta_a.^2).*alph_ms.*...
                  alph_ar.^(-kc);
              
bet_xs1 = bet_ms-0.5.*bet_ar;     
alph_xs1 = (2) .*g.*(rhoa./eta_a.^2).*alph_ms.*...
                  alph_ar.^(-0.5);
              
Xnew0 =(2.*alph_ms.*g.*rhoa./(sig.*eta_a.^2)).*D.^(bet_xs);
Xnew1 = alph_xs1.*D.^(bet_xs1);

[bm_new,am_new] = fallspeed_params(0,Xnew0);
[bm_new2,am_new2] = fallspeed_params(2,Xnew0);
[bm_new3,am_new3] = fallspeed_params(3,Xnew0);
[bm_new1,am_new1] = fallspeed_params(1,Xnew1);

%alph_vt for area ratio param.
%alph_v_sph = 0.5.*(eta_a./rhoa).*am_sph.*alph_xs.^bm_sph;
%alph_v_0 = 0.5.*(eta_a./rhoa).*am_new.*alph_xs.^bm_new;
alph_v_0 = am_new.*nua.^(1-2.*bm_new).*(2.*alph_ms.*g./(rhoa.*sig)).^bm_new;
bet_v_0 = bm_new.*bet_xs-1;

alph_v_2 = am_new2.*nua.^(1-2.*bm_new2).*(2.*alph_ms.*g./(rhoa.*sig)).^bm_new2;
bet_v_2 = bm_new2.*bet_xs-1;

alph_v_3 = am_new3.*nua.^(1-2.*bm_new3).*(2.*alph_ms.*g./(rhoa.*sig)).^bm_new3;
bet_v_3 = bm_new3.*bet_xs-1;

alph_v_1 = am_new1.*nua.^(1-2.*bm_new2).*(2.*alph_ms.*g./(rhoa.*sig)).^bm_new1;
bet_v_1 = bm_new1.*bet_xs1-1;


vt0 = alph_v_0.*D.^bet_v_0;
vt2 = alph_v_2.*D.^bet_v_2;
vt3 = alph_v_3.*D.^bet_v_3;
vt1 = alph_v_1.*D.^bet_v_1;

figure; 
ax = tight_subplot(3,2,[0 0.1],[0.1 0.02],[0.1 0.02]);
%subplot(3,2,2)
plot(ax(2),D,vt2./100,'k','linewidth',2.0)
hold(ax(2));
plot(ax(2),D,vt0./100,'--k','linewidth',2.0);
plot(ax(2),D,vt1./100,':k','linewidth',2.0)
plot(ax(2),D,vt3./100,'-.k','linewidth',2.0);
legend(ax(2),'M96/KC02','MH05','HW10','KC05','interpreter','latex','fontsize',16,'location','NorthWest')
ylim(ax(2),[0 2])
set(ax(2),'xscale','log','xticklabel','','ticklabelinterpreter','latex','FontSize',fsize,...
    'xgrid','on','ygrid','on','xminortick','on','yminortick','on','linewidth',1.4,...
    'ticklength',[0.025 0.05]);

%subplot(3,2,4)
plot(ax(4),D,bet_v_2,'k','linewidth',2.0)
hold(ax(4))
plot(ax(4),D,bet_v_0,'--k','linewidth',2.0);
plot(ax(4),D,bet_v_1,':k','linewidth',2.0)
plot(ax(4),D,bet_v_3,'-.k','linewidth',2.0);
set(ax(4),'xscale','log','xticklabel','','ticklabelinterpreter','latex','FontSize',fsize,...
    'xgrid','on','ygrid','on','xminortick','on','yminortick','on','linewidth',1.4,...
    'ticklength',[0.025 0.05]);

%subplot(3,2,6)
plot(ax(6),D,alph_v_2,'k','linewidth',2.0)
hold(ax(6))
plot(ax(6),D,alph_v_0,'--k','linewidth',2.0);
plot(ax(6),D,alph_v_1,':k','linewidth',2.0)
plot(ax(6),D,alph_v_3,'-.k','linewidth',2.0);
set(ax(6),'xscale','log','ticklabelinterpreter','latex','FontSize',fsize,...
    'xgrid','on','ygrid','on','xminortick','on','yminortick','on','linewidth',1.4,...
    'ticklength',[0.025 0.05]);


Re0 = am_0.*X.^bm_0;
Re1 = am_1.*X.^bm_1;
Re2 = am_2.*X.^bm_2;
Re3 = am_3.*X.^bm_3;

%subplot(3,2,1)
plot(ax(1),X,Re2,'k','linewidth',2.0);
hold(ax(1))
plot(ax(1),X,Re0,'--k','linewidth',2.0);
plot(ax(1),X,Re1,':k','linewidth',2.0);
plot(ax(1),X,Re3,'-.k','linewidth',2.0);
%title(ax(1),'X-Re')
xlim(ax(1),10.^([1 8]));
set(ax(1),'xscale','log','yscale','log','xticklabel','','ticklabelinterpreter','latex','FontSize',fsize,...
    'xgrid','on','ygrid','on','xminortick','on','yminortick','on','linewidth',1.4,...
    'ticklength',[0.025 0.05]);

%subplot(3,2,3)
plot(ax(3),X,bm_2,'k','linewidth',2.0);
hold(ax(3))
plot(ax(3),X,bm_0,'--k','linewidth',2.0);
plot(ax(3),X,bm_1,':k','linewidth',2.0);
plot(ax(3),X,bm_3,'-.k','linewidth',2.0);
xlim(ax(3),10.^([1 8]));
%title(ax(3),'X-bm')
set(ax(3),'xscale','log','xticklabel','','ticklabelinterpreter','latex','FontSize',fsize,...
    'xgrid','on','ygrid','on','xminortick','on','yminortick','on','linewidth',1.4,...
    'ticklength',[0.025 0.05]);

%subplot(3,2,5)
plot(ax(5),X,am_2,'k','linewidth',2.0);
hold(ax(5))
plot(ax(5),X,am_0,'--k','linewidth',2.0);
plot(ax(5),X,am_1,':k','linewidth',2.0);
plot(ax(5),X,am_3,'-.k','linewidth',2.0);
%title(ax(5),'X-am');
xlim(ax(5),10.^([1 8]));

set(ax(5),'xscale','log','yscale','log','ticklabelinterpreter','latex','FontSize',fsize,...
    'xgrid','on','ygrid','on','xminortick','on','yminortick','on','linewidth',1.4,...
    'ticklength',[0.025 0.05]);




function [bm,am] = fallspeed_params(Heyms_switch,Xbar)
% Aggregate fallspeed parameters 

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

%bm = ...
%(C1.*sqrt(Xbar)./...
%(2.*(sqrt(1+C1.*sqrt(Xbar)) - 1).*sqrt(1+C1.*sqrt(Xbar)))) -...
%((ao.*bo.*Xbar.^bo)./(C2.*(sqrt(1+C1.*sqrt(Xbar))-1).^2)) + ...
%dbm;

bml = (C1.*sqrt(Xbar)./...
(2.*(sqrt(1+C1.*sqrt(Xbar)) - 1).*sqrt(1+C1.*sqrt(Xbar)))) -...
((ao.*bo.*Xbar.^bo)./(C2.*(sqrt(1+C1.*sqrt(Xbar))-1).^2));

bm = bml+dbm;

am = ...
Psi_tX.*(C2.*((sqrt(1+C1.*sqrt(Xbar))-1).^2) - ao.*Xbar.^bo)./(Xbar.^(bml));


end