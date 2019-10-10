function [lik, latents] = R4(param,data)
    
    % Rationally inattentive discounting with endogenized temperature.
    
    beta = param(1);           % magnitude effect parameter
    k0 = param(2);             % baseline discounting
    alpha0 = param(3);         % baseline inverse temperature
    
    if length(param) > 3
        lapse = param(4);
    else
        lapse = 0;
    end
    
    r = [data.X2 data.X1];  % rewards
    t = [data.T2 data.T1];  % delays
    s2_e = data.s2_u./(abs(r).*beta);
    k = max(0,k0 + s2_e./data.s2_u);
    D = 1./(1 + (k.*t));
    V = r.*D;    % discounted value
    alpha = max(1e-6,alpha0 + 1./sqrt(sum(((D.^2).*s2_e),2)));
    P = (1-lapse)*normcdf(alpha.*(V(:,1)-V(:,2))) + lapse*0.5; % choice probability
    lik = safelog(P')*data.LL + safelog(1-P')*(1-data.LL);  % log likelihood
    
    if nargout > 1
        latents.P = P;
        latents.V = V;
        latents.alpha = alpha;
        latents.k = k;
        latents.s2_e = s2_e;
    end