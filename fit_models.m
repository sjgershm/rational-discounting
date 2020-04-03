function [results, bms_results] = fit_models(data,models,results)
    
    % Fit models to data.
    %
    % USAGE: [results, bms_results] = fit_models(data,models,results)
    
    likfuns = {'QH' 'H0' 'H1' 'H2' 'H3' 'R1' 'R2'};
    if nargin < 2; models = 1:length(likfuns); end
    
    mb = [-5 5];
    cb = [-5 5];
    
    for i = models
        
        clear param
        switch likfuns{i}
            
            case 'QH'
                param(1) = struct('name','alpha','lb',1e-6,'ub',50);
                param(2) = struct('name','delta','lb',1e-6,'ub',1);
                param(3) = struct('name','beta','lb',1e-6,'ub',1);
                param(4) = struct('name','lapse','lb',0,'ub',0.99);
            
            case 'H0'
                param(1) = struct('name','alpha','lb',1e-6,'ub',50);
                param(2) = struct('name','k','lb',1e-6,'ub',5);
                param(3) = struct('name','lapse','lb',0,'ub',0.99);
            
            case 'H1'
                param(1) = struct('name','alpha','lb',1e-6,'ub',50);
                param(2) = struct('name','m_k','lb',-5,'ub',5);
                param(3) = struct('name','c_k','lb',-5,'ub',5);
                param(4) = struct('name','lapse','lb',0,'ub',0.99);
                
            case 'H2'
                param(1) = struct('name','m_alpha','lb',mb(1),'ub',mb(2));
                param(2) = struct('name','c_alpha','lb',cb(1),'ub',cb(2));
                param(3) = struct('name','m_k','lb',mb(1),'ub',mb(2));
                param(4) = struct('name','c_k','lb',cb(1),'ub',cb(2));
                param(5) = struct('name','lapse','lb',0,'ub',0.99);
                
            case 'H3'
                param(1) = struct('name','m_alpha','lb',mb(1),'ub',mb(2));
                param(2) = struct('name','m_k','lb',mb(1),'ub',mb(2));
                param(3) = struct('name','lapse','lb',0,'ub',0.99);
                
            case 'R1'
                param(1) = struct('name','alpha','lb',1e-6,'ub',50);
                param(2) = struct('name','beta','lb',1e-6,'ub',100);
                param(3) = struct('name','lapse','lb',0,'ub',0.99);
                
            case 'R2'
                param(1) = struct('name','beta','lb',1e-6,'ub',100);
                param(2) = struct('name','lapse','lb',0,'ub',0.99);
                
        end
        
        disp(['... ',likfuns{i}]);
        results(i) = mfit_optimize(str2func(likfuns{i}),param,data);
    end
    
    if nargout > 1
        bms_results = mfit_bms(results,1);
    end