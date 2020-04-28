function out = Fox_H(bm, Bm, an, An, ap, Ap, bq, Bq, z)
% This matlab function is given in:
% "Performance Analysis of M-QAM Multihop Relaying over mmWave Weibull
%  Fading Channels", Soulimani et al. 2016
    
abstol = 1e-12;
reltol = 1e-8;

    
%% Integrand definition

F = @(s,z) (GammaProd(bm,Bm,s) .* GammaProd(1-an,-An,s) .* z.^(-s)) ./ ...
         (GammaProd(1-bq,-Bq,s) .* GammaProd(ap,Ap,s));
     
%% Contour Preparation:
epsilon = 1;

Sups = min((1-an)./An); Infs = max(-bm./Bm);
if(isempty(Sups) && isempty(Infs))
    WPx=1;
elseif(isempty(Sups) && ~isempty(Infs))
    WPx = Infs + epsilon;
elseif(~isempty(Sups) && isempty(Infs))
    WPx = Sups-epsilon;
else 
    WPx = (Sups + Infs)/2; % s between Sups and Infs
end

%% Integration
infity = 100;
    out = (1./(2i.*pi)).*integral(@(s) F(s,z),WPx-1i.*infity,WPx+1i.*infity,'arrayvalued',true,...
                'abstol',abstol,'reltol',reltol);

out = real(out);

return





end