function [lik, latents] = H2(param,data)
    
    % Hyperbolic discounting with magnitude-dependent discount rate and
    % temperature.
    
    m_alpha = param(1);
    c_alpha = param(2);
    m_k = param(3);
    c_k = param(4);
    
    if length(param) > 4
        lapse = param(5);
    else
        lapse = 0;
    end
    
    r = [data.X2 data.X1];  % rewards
    t = [data.T2 data.T1];  % delays
    k = exp(c_k - m_k.*log(abs(r)));
    alpha = exp(c_alpha + m_alpha.*log(abs(r(:,1))+abs(r(:,2))));
    V = r./(1 + (k.*t));    % discounted value
    P = (1-lapse)*normcdf(alpha.*(V(:,1)-V(:,2))) + lapse*0.5; % choice probability
    lik = safelog(P')*data.LL + safelog(1-P')*(1-data.LL);  % log likelihood
    
    if nargout > 1
        latents.P = P;
        latents.Vdiff = V(:,1)-V(:,2);
    end