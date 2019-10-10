function [m,se,X] = quantile_stats(x,y,q)
    
    if isscalar(q)
        q = quantile(y,q);
    end
    
    for i = 1:length(q)-1
        ix = y>q(i) & y<=q(i+1);
        X{i} = x(ix);
        m(i) = nanmean(x(ix));
        se(i) = nanstd(x(ix))./sqrt(sum(~isnan(x(ix))));
    end