function [results, bms_results] = fit_models(data,models,results)
    
    % Fit models to data.
    %
    % USAGE: [results, bms_results] = fit_models(data,models,results)
    
    likfuns = {'M1' 'M2' 'M3' 'M4' 'M5'};
    if nargin < 2; models = 1:length(likfuns); end
    
    for i = models
        
        clear param
        switch likfuns{i}
            
            case 'M1'
                param(1) = struct('name','alpha','lb',1e-6,'ub',50);
                param(2) = struct('name','m_k','lb',-5,'ub',5);
                param(3) = struct('name','c_k','lb',-5,'ub',5);
                param(4) = struct('name','lapse','lb',0,'ub',0.99);
                
            case 'M2'
                param(1) = struct('name','m_alpha','lb',-5,'ub',5);
                param(2) = struct('name','c_alpha','lb',-5,'ub',5);
                param(3) = struct('name','m_k','lb',-5,'ub',5);
                param(4) = struct('name','c_k','lb',-5,'ub',5);
                param(5) = struct('name','lapse','lb',0,'ub',0.99);
                
            case 'M3'
                param(1) = struct('name','m_alpha','lb',-5,'ub',5);
                param(2) = struct('name','m_k','lb',-5,'ub',5);
                param(3) = struct('name','lapse','lb',0,'ub',0.99);
                
            case 'M4'
                param(1) = struct('name','alpha','lb',1e-6,'ub',50);
                param(2) = struct('name','beta','lb',1e-6,'ub',100);
                param(3) = struct('name','lapse','lb',0,'ub',0.99);
                
            case 'M5'
                param(1) = struct('name','beta','lb',1e-6,'ub',100);
                param(2) = struct('name','lapse','lb',0,'ub',0.99);
                
        end
        
        disp(['... ',likfuns{i}]);
        results(i) = mfit_optimize(str2func(likfuns{i}),param,data);
    end
    
    if nargout > 1
        bms_results = mfit_bms(results);
    end