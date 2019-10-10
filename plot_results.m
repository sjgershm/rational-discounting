function plot_results(fig,results,data,varargin)
    
    rng(2); % set random seed for reproducibility
    
    switch fig
        
        case 'fig1'
            
            d = linspace(0,20,100);
            figure;
            v = 0.3*d;
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
            legend({'Likelihood' 'Prior' 'Posterior'},'FontSize',25,'Location','South','Box','Off')
            
        case 'fig2'
            
            if nargin < 2; load results/results_hardisty; end
            if nargin < 3; data = load_data('hardisty'); end
            
            cond = [data.Condition];
            
            r = [10 100 1000 10000];
            for i = 1:4
                ix = cond==i;
                beta = mean(results(7).x(ix,1));
                k(:,i) = [1./(beta*r(i)); -1./(beta*r(i))];
            end
            
            bar(k);
            set(gca,'YLim',[min(k(:))-std(k(:))./2 max(k(:))+std(k(:))/2],'FontSize',25,'XTickLabel',{'Gain' 'Loss'});
            xlabel('Reward magnitude','FontSize',25);
            ylabel('Discount factor (k)','FontSize',25);
            legend({'$10' '$100' '$1000' '$10000'},'FontSize',25,'Box','Off');
            
        case 'fig3'
            
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
            
            [se,m] = wse(LL);
            x = q(1:end-1) + diff(q)/2;
            clr = [0 0 0; 0.5 0.5 0.5];
            mm = squeeze(nanmean(LLm));
            for i = 1:2
                errorbar(x',m(:,i),se(:,i),'o','MarkerFaceColor',clr(i,:),'MarkerSize',10,'Color',clr(i,:));
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
            se = std(results.x(:,1:2))./sqrt(size(results.x,1));
            
            barerrorbar(m',se'); colormap summer
            set(gca,'Fontsize',25,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'m_\alpha' 'm_k'});
            ylabel('Parameter estimate','FontSize',25);
            [~,p,~,stat] = ttest(results.x(:,1:2));
            disp(['t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            
        case 'chavez_corr'
            
            N = 7;
            [m,se] = quantile_stats(results(5).x(:,1),results(5).x(:,2),N);
            p = linspace(0,1,N);
            p = p(1:end-1) + diff(p)/2;
            errorbar(p,m,se,'-ok','LineWidth',4,'MarkerFaceColor','k','MarkerSize',12);
            xlabel('m_k quantile','FontSize',25);
            ylabel('m_\alpha','FontSize',25);
            set(gca,'FontSize',25,'XLim',[0 1]);
            
        case 'fig4'
            
            data = load_data('bhui');
            load results/results_bhui
            
            results = results(5);
            
            for i = 1:2
                if i==1; c=[data.Condition]==0; else c=[data.Condition]==1; end
                m(i,:) = mean(results.x(c,1:2));
                se(i,:) = std(results.x(c,1:2))./sqrt(sum(c));
            end
            
            c = [data.Condition]==0;
            [~,p,~,stat] = ttest2(results.x(c,1),results.x(~c,1));
            disp(['m_alpha (low vs. high): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            [~,p,~,stat] = ttest2(results.x(c,3),results.x(~c,3));
            disp(['m_k (low vs. high): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            [~,p,~,stat] = ttest(results.x(:,1));
            disp(['m_alpha: t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            [~,p,~,stat] = ttest(results.x(:,3));
            disp(['m_k: t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            
            barerrorbar(m',se'); colormap summer;
            set(gca,'Fontsize',25,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'m_\alpha' 'm_k'});
            ylabel('Parameter estimate','FontSize',25);
            
            legend({'Low variance' 'High variance'},'FontSize',25,'Box','Off');
            
    end