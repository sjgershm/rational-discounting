function [lik, latents] = M3(param,data)
    
    m_alpha = param(1);
    m_k = param(2);
    
    if length(param) > 2
        lapse = param(3);
    else
        lapse = 0;
    end
    
    r = [data.X2 data.X1];  % rewards
    t = [data.T2 data.T1];  % delays
    k = exp(m_k.*log(r));
    alpha = exp(m_alpha.*log(r(:,1)+r(:,2)));
    V = r./(1 + (k.*t));    % discounted value
    P = (1-lapse)*normcdf(alpha.*(V(:,1)-V(:,2))) + lapse*0.5; % choice probability
    lik = safelog(P')*data.LL + safelog(1-P')*(1-data.LL);  % log likelihood
    
    if nargout > 1
        latents.P = P;
        latents.Vdiff = V(:,1)-V(:,2);
    end