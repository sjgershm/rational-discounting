function simdata = simulate_data(data,results)
    
    % Simulate the R2 model on Bhui data.
    
    for s = 1:length(data)
        simdata(s) = R2_sim(results(7).x(s,:),data(s));
    end