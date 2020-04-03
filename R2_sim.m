function data = R2_sim(param,data)
    
    % Rationally inattentive discounting with endogenized temperature.
    
    beta = param(1);           % magnitude effect parameter
    
    if length(param) > 1
        lapse = param(2);
    else
        lapse = 0;
    end
    
    r = [data.X2 data.X1];  % rewards
    t = [data.T2 data.T1];  % delays
    s2_e = data.s2_u./(abs(r).*beta);
    k = s2_e./data.s2_u;
    D = 1./(1 + (k.*t));
    V = r.*D;    % discounted value
    alpha = 1./sqrt(sum(((D.^2).*s2_e),2));
    P = (1-lapse)*normcdf(alpha.*(V(:,1)-V(:,2))) + lapse*0.5; % choice probability
    
    data.LL = double(rand(size(P)) < P);