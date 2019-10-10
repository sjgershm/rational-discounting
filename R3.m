function [lik, latents] = R3(param,data)
    
    % Rationally inattentive discounting.
    
    alpha = param(1);       % inverse temperature
    beta = param(2);        % magnitude effect parameter
    k0 = param(3);          % baseline disounting
    
    if length(param) > 3
        lapse = param(4);
    else
        lapse = 0;
    end
    
    r = [data.X2 data.X1];  % rewards
    t = [data.T2 data.T1];  % delays
    s2_e = data.s2_u./(abs(r).*beta);
    k = max(0,k0 + s2_e./data.s2_u);
    V = r./(1 + (k.*t));
    P = (1-lapse)*normcdf(alpha*(V(:,1)-V(:,2))) + lapse*0.5; % choice probability
    lik = safelog(P')*data.LL + safelog(1-P')*(1-data.LL);  % log likelihood
    
    if nargout > 1
        latents.P = P;
        latents.Vdiff = V(:,1)-V(:,2);
    end