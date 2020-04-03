function [lik, latents] = H2_lik(param,data)
    
    % Hyperbolic discounting.
    
    alpha1 = param(1);       % inverse temperature (low reward)
    alpha2 = param(2);       % inverse temperature (high reward)
    k1 = param(3);           % discount parameter (low reward)
    k2 = param(4);           % discount parameter (high reward)
    
    if length(param) > 4
        lapse = param(5);
    else
        lapse = 0;
    end
    
    ix = data.X2 > median(data.X2);
    k = zeros(size(ix)) + k1;
    k(ix) = k2;
    alpha = zeros(size(ix)) + alpha1;
    alpha(ix) = alpha2;
    
    r = [data.X2 data.X1];  % rewards
    t = [data.T2 data.T1];  % delays
    V = r./(1 + (k.*t));    % discounted value
    P = (1-lapse)*normcdf(alpha.*(V(:,1)-V(:,2))) + lapse*0.5; % choice probability
    lik = safelog(P')*data.LL + safelog(1-P')*(1-data.LL);  % log likelihood
    
    if nargout > 1
        latents.P = P;
        latents.V = V;
    end