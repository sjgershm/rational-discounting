function [lik, latents] = M1(param,data)
    
    alpha = param(1);       % inverse temperature
    m_k = param(2);
    c_k = param(3);
    
    if length(param) > 3
        lapse = param(4);
    else
        lapse = 0;
    end
    
    r = [data.X2 data.X1];  % rewards
    t = [data.T2 data.T1];  % delays
    k = exp(c_k + m_k.*log(r));
    V = r./(1 + (k.*t));    % discounted value
    P = (1-lapse)*normcdf(alpha.*(V(:,1)-V(:,2))) + lapse*0.5; % choice probability
    lik = safelog(P')*data.LL + safelog(1-P')*(1-data.LL);  % log likelihood
    
    if nargout > 1
        latents.P = P;
        latents.V = V;
    end