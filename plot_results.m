function plot_results(fig,results,data,varargin)
    
    rng(2); % set random seed for reproducibility
    
    switch fig
        
        case 'fig1'
            
            d = linspace(0,20,100);
            figure;
            s = [0.3 0.03];
            T = {'A' 'B'};
            for i = 1:2
                subplot(2,2,i);
                v = s(i)*d;
                shadedErrorBar(d,ones(size(d)),sqrt(v));
                hold on;
                x = linspace(-4,4,100);
                y = normpdf(x,0,1);
                plot(y*5,x,'LineWidth',4);
                w = 1./(1+v);
                shadedErrorBar(d,w,sqrt(1./((1./v)+1)),'LineProps','-b');
                set(gca,'FontSize',25);
                xlabel('Delay','FontSize',25);
                ylabel('Value','FontSize',25);
                if i==1
                    legend({'Likelihood' 'Prior' 'Posterior'},'FontSize',25,'Location','South','Box','Off');
                end
                mytitle(T{i},'Left','FontSize',25,'FontWeight','Bold')
            end
            
            subplot(2,2,3)
            b = linspace(0.5,5,100);
            r = [1 2 3];
            C = linspecer(length(r));
            for i = 1:length(r)
                plot(b,1./(b.*r(i)),'Color',C(i,:),'LineWidth',4);
                hold on;
            end
            set(gca,'FontSize',25,'XLim',[min(b) max(b)])
            legend({'Low reward' 'Medium reward' 'High reward'},'FontSize',25,'Box','Off')
            xlabel('Sensitivity parameter (\beta)','FontSize',25);
            ylabel('Discount factor (k)','FontSize',25);
            mytitle('C','Left','FontSize',25,'FontWeight','Bold')
            
            subplot(2,2,4)
            b = linspace(0.5,5,100);
            t = 0.1;
            C = linspecer(length(r));
            for i = 1:length(r)
                alpha = sqrt(b*r(i)/t).*(1 + t./(b.*r(i)));
                plot(b,1./alpha,'Color',C(i,:),'LineWidth',4);
                hold on;
            end
            set(gca,'FontSize',25,'XLim',[min(b) max(b)])
            xlabel('Sensitivity parameter (\beta)','FontSize',25);
            ylabel('Choice stochasticity (1/\alpha)','FontSize',25);
            mytitle('D','Left','FontSize',25,'FontWeight','Bold')
            
            set(gcf,'Position',[200 200 1000 800])
            
        case 'fig2'
            
            rng(2);
            load ballard_data_exp3
            M = [20 50 100 200 2000];
            zs = abs(zscore(data.choice));
            K = []; X = []; Z = []; xx = [];
            for j = 1:2
                for i = 1:5
                    ix = data.magnitude==i & data.manipulation==j & zs<4;
                    c = data.choice(ix);
                    k = ((c./M(i))-1)/30;
                    K = [K; k];
                    n = length(k);
                    X = [X; zeros(n,1)-j+1 zeros(n,1)+M(i)];
                    m(i,j) = mean(k);
                    ci = bootci(2000,@mean,k);
                    se(i,j) = diff(ci)/2;
                    z(i,j) = std(k);
                    Z = [Z; z(i,j)];
                    xx = [xx; 1-j M(i)];
                    ci = bootci(2000,@std,k);
                    se_z(i,j) = diff(ci)/2;
                end
            end
            
            [b,~,stat] = glmfit([X X(:,1).*X(:,2)],K);
            b
            stat.p
            
            [b,~,stat] = glmfit([xx xx(:,1).*xx(:,2)],Z);
            b
            stat.p
            
            figure;
            subplot(1,2,1);
            errorbar(m,se,'-o','LineWidth',4,'MarkerFaceColor','w','MarkerSize',10)
            set(gca,'FontSize',25,'XTickLabel',M,'XLim',[0.5 5.5]);
            xlabel('Smaller sooner magnitude','FontSize',25);
            ylabel('Discount factor (k)','FontSize',25);
            legend({'No justification' 'Justification'},'FontSize',25,'Box','Off');
            mytitle('A','Left','FontSize',25,'FontWeight','Bold')
            
            subplot(1,2,2);
            errorbar(z,se_z,'-o','LineWidth',4,'MarkerFaceColor','w','MarkerSize',10)
            set(gca,'FontSize',25,'XTickLabel',M,'XLim',[0.5 5.5]);
            xlabel('Smaller sooner magnitude','FontSize',25);
            ylabel('Standard deviation','FontSize',25);
            mytitle('B','Left','FontSize',25,'FontWeight','Bold')
            
            set(gcf,'Position',[200 200 950 400])
            
        case 'fig3'
            
            if nargin < 3; data = load_data('chavez'); end
            if nargin < 2; load results/results_chavez; end
            
            subplot(2,2,1)
            plot_results('chavez_parameters',results,data)
            mytitle('A','Left','FontSize',25,'FontWeight','Bold')
            subplot(2,2,2)
            plot_results('chavez_corr',results,data)
            mytitle('B','Left','FontSize',25,'FontWeight','Bold')
            subplot(2,2,3)
            plot_results('psychometric',results,data,2)
            mytitle('C                  Hyperbolic','Left','FontSize',25,'FontWeight','Bold')
            subplot(2,2,4)
            plot_results('psychometric',results,data,7)
            mytitle('D                  Rational','Left','FontSize',25,'FontWeight','Bold')
            legend off;
            
            set(gcf,'Position',[200 200 1200 900])
            
            
        case 'psychometric'
            
            if nargin < 4; model = 7; else model = varargin{1}; end
            if nargin < 3; data = load_data('chavez'); end
            if nargin < 2; load results/results_chavez; end
            
            q = quantile(data(1).T2,7);
            r = data(1).X1 + data(1).X2;
            
            for s = 1:length(data)
                ix = r < median(r);
                LL(s,:,1) = quantile_stats(data(s).LL(ix),data(s).T2(ix),q);
                LLm(s,:,1) = quantile_stats(results(model).latents(s).P(ix),data(s).T2(ix),q);
                ix = r >= median(r);
                LL(s,:,2) = quantile_stats(data(s).LL(ix),data(s).T2(ix),q);
                LLm(s,:,2) = quantile_stats(results(model).latents(s).P(ix),data(s).T2(ix),q);
            end
            
            %[se,m] = wse(LL);
            m = squeeze(nanmean(LL));
            x = q(1:end-1) + diff(q)/2;
            clr = [0 0 0; 0.5 0.5 0.5];
            mm = squeeze(nanmean(LLm));
            for i = 1:2
                %errorbar(x',m(:,i),se(:,i),'o','MarkerFaceColor',clr(i,:),'MarkerSize',10,'Color',clr(i,:));
                plot(x',m(:,i),'o','MarkerFaceColor',clr(i,:),'MarkerSize',12,'Color',clr(i,:));
                hold on;
                plot(x',mm(:,i),'-','LineWidth',4,'Color',clr(i,:));
            end
            legend({'Data: low mag' 'Model: low mag' 'Data: high mag' 'Model: high mag'},'FontSize',25,'Box','Off');
            set(gca,'FontSize',25,'YLim',[0 1]);
            xlabel('Delay (days)','FontSize',25);
            ylabel('P(larger later)','FontSize',25);
            
        case 'chavez_parameters'
            
            if nargin < 3; data = load_data('chavez'); end
            if nargin < 2; load results/results_chavez; end
            
            results = results(5);
            
            m = mean(results.x(:,1:2));
            ci = bootci(2000,@mean,results.x(:,1:2));
            err= diff(ci)/2;
            %se = std(results.x(:,1:2))./sqrt(size(results.x,1));
            
            barerrorbar(m',err'); colormap summer
            set(gca,'Fontsize',25,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'m_\alpha' 'm_k'});
            ylabel('Parameter estimate','FontSize',25);
            [~,p,~,stat] = ttest(results.x(:,1:2));
            disp(['t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            
        case 'chavez_corr'
            
            N = 7;
            [m,err] = quantile_stats(results(5).x(:,1),results(5).x(:,2),N);
            p = linspace(0,1,N);
            p = p(1:end-1) + diff(p)/2;
            errorbar(p,m,err,'-ok','LineWidth',4,'MarkerFaceColor','k','MarkerSize',12);
            xlabel('m_k quantile','FontSize',25);
            ylabel('m_\alpha','FontSize',25);
            set(gca,'FontSize',25,'XLim',[0 1]);
            
        case 'fig4'
            
            rng(1);
            
            if nargin < 2; load results/results_bhui; end
            if nargin < 3; data = load_data('bhui'); end
            
            %b = results(7).x(:,1);
            %for s=1:length(data); z(s) = mean(b(s).*data(s).X2>data(s).T2); end
            %ix = results(4).x(:,3)'>0;
            %ix = z >=0;
            
            results = results(4);
            
            X = [];
            for i = 1:2
                if i==1; c=[data.Condition]==0; else c=[data.Condition]==1; end
                %c = c & ix;
                m(i,:) = mean(results.x(c,[1 3]));
                X = [X; results.x(c,[1 3])];
                %err(i,:) = std(results.x(c,[1 3]))./sqrt(sum(c));
            end
            
            ci = bootci(2000,@mean,X); err = repmat(diff(ci)/2,2,1);    % pooled variance
            
            c = [data.Condition]==0;
            %[~,p,~,stat] = ttest2(results.x(c & ix,1),results.x(~c & ix,1));
            [~,p,~,stat] = ttest2(results.x(c,1),results.x(~c,1));
            disp(['m_alpha (low vs. high): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            %[~,p,~,stat] = ttest2(results.x(c & ix,3),results.x(~c & ix,3));
            [~,p,~,stat] = ttest2(results.x(c,3),results.x(~c,3));
            disp(['m_k (low vs. high): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            %[~,p,~,stat] = ttest(results.x(ix,1));
            [~,p,~,stat] = ttest(results.x(:,1));
            disp(['m_alpha: t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            %[~,p,~,stat] = ttest(results.x(ix,3));
            [~,p,~,stat] = ttest(results.x(:,3));
            disp(['m_k: t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            
            barerrorbar(m',err'); colormap summer;
            set(gca,'Fontsize',25,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'m_\alpha' 'm_k'});
            ylabel('Parameter estimate','FontSize',25);
            legend({'Low variance' 'High variance'},'FontSize',25,'Box','Off');
            
    end