function [lik, latents] = R3_lik(param,data)
    
    % Rationally inattentive discounting with endogenized temperature.
    
    k0 = param(1);             % baseline discounting
    b_k = param(2);             % discounting magnitude effect parameter
    alpha0 = param(3);         % baseline inverse temperature
    b_alpha = param(4);         % choice stochasticity magnitude effect parameter
    %s2_u = param(5);            % reward variance
    
    if length(param) > 4
        lapse = param(5);
    else
        lapse = 0;
    end
    
    r = [data.X2 data.X1];  % rewards
    t = [data.T2 data.T1];  % delays
    k = k0 + 1./(abs(r).*b_k);
    D = 1./(1 + (k.*t));
    V = r.*D;    % discounted value
    
    s2_u = 1;
    s2_e = s2_u./(abs(r).*b_alpha);
    k = s2_e./s2_u;
    D_alpha = 1./(1 + (k.*t));
    alpha = alpha0 + 1./sqrt(sum(((D_alpha.^2).*s2_e.*t),2));
    P = (1-lapse)*normcdf(alpha.*(V(:,1)-V(:,2))) + lapse*0.5; % choice probability
    lik = safelog(P')*data.LL + safelog(1-P')*(1-data.LL);  % log likelihood
    
    if nargout > 1
        latents.P = P;
        latents.V = V;
        latents.alpha = alpha;
        latents.k = k;
        latents.s2_e = s2_e;
    end