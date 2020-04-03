function [m,err,X] = quantile_stats(x,y,q)
    
    rng(1);
    
    if isscalar(q)
        q = quantile(y,q);
    end
    
    for i = 1:length(q)-1
        ix = y>q(i) & y<=q(i+1);
        X{i} = x(ix);
        m(i) = nanmean(x(ix));
        if nargout > 1
            ci = bootci(2000,@nanmean,x(ix));
            err(i) = diff(ci)/2;
        end
        %se(i) = nanstd(x(ix))./sqrt(sum(~isnan(x(ix))));
    end