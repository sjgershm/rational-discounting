function [ci,bootstat] = ballard_regression
    
    load ballard_data_exp3
    
    for i = 1:length(data.condition)
        d(i).condition = data.condition(i);
        d(i).manipulation = data.manipulation(i);
        d(i).magnitude = data.magnitude(i);
        d(i).choice = data.choice(i);
    end
    
    zs = abs(zscore([d.choice]));
    d(zs>4) = [];
    
    [ci,bootstat] = bootci(1000,@regressfun,d);
    
end

function b = regressfun(d)
    
    M = [20 50 100 200 2000];
    xx = []; MM = []; Z = [];
    v = [1 0];
    
    for j = 1:2
        for i = 1:5
            ix = [d.magnitude]==i & [d.manipulation]==j;
            c = [d(ix).choice];
            k = ((c./M(i))-1)/30;
            Z = [Z; std(k)];
            MM = [MM; mean(k)];
            %xx = [xx; -v(j) -i];
            xx = [xx; v(j) i];
        end
    end
    
    x = [xx xx(:,1).*xx(:,2)];
    b1 = glmfit(x,MM);
    b2 = glmfit(x,Z);
    b = [b1' b2'];
    
end