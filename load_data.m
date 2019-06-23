function data = load_data(dataset)
    
    % Load data from csv file.
    %
    % USAGE: data = load_data(dataset)
    
    switch dataset
            
        case 'chavez'
            
            D = csvread('data/chavez_data.csv',1);
            
            X1 = [54 55 19 31 14 47 15 25 78 40 11 67 34 27 69 49 80 24 33 28 34 25 41 54 54 22 20]';
            T1 = zeros(27,1);
            X2 = [55 75 25 85 25 50 35 60 80 55 30 75 35 50 85 60 85 35 80 30 50 30 75 60 80 25 55]';
            T2 = [117 61 53 7 19 160 13 14 162 62 7 119 186 21 91 89 157 29 14 179 30 80 20 111 30 136 7]';
            
            for s = 1:size(D,1)
                data(s).LL = D(s,3:end)';
                data(s).X1 = X1;
                data(s).X2 = X2;
                data(s).T1 = T1;
                data(s).T2 = T2;
                data(s).N = length(X1);
                data(s).C = 2;
                data(s).s2_u = var([data(s).X1; data(s).X2]);
            end
            
        case 'bhui'
            
            % Condition = 0: low variance
            % Condition = 1: high variance
            
            D = readtable('data/bhui_data.csv');
            S = unique(D.Subject);
            F = fieldnames(D); F = F(1:7);
            
            for s = 1:length(S)
                ix = D.Subject==S(s);
                for f = 1:length(F)
                    data(s).(F{f}) = D.(F{f})(ix,:);
                end
                data(s).Condition = data(s).Condition(1);
                data(s).LL = D.LaterOptionChosen(ix,:);
                bad = isnan(data(s).LL);
                data(s).X1(bad) = []; data(s).X2(bad) = [];
                data(s).T1(bad) = []; data(s).T2(bad) = [];
                data(s).LL(bad) = [];
                data(s).N = length(data(s).LL);
                data(s).C = 2;
                data(s).s2_u = var([data(s).X1; data(s).X2]);
            end
    end