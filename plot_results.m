function plot_results(fig,results,data)
    
    rng(2); % set random seed for reproducibility
    
    switch fig
        
        case 'bhui_parameters'
            
            if nargin < 3; data = load_data('bhui'); end
            if nargin < 2; load results/results_bhui; end
            
            results(3).x(:,2) = -results(3).x(:,2); % flip sign on discounting magnitude effect
            
            for i = 1:2
                if i==1; c=[data.Condition]==0; else c=[data.Condition]==1; end
                m(i,:) = mean(results(3).x(c,1:2));
                se(i,:) = std(results(3).x(c,1:2))./sqrt(sum(c));
            end
            
            c = [data.Condition]==0;
            [~,p,~,stat] = ttest2(results(3).x(c,1),results(3).x(~c,1));
            disp(['m_alpha (low vs. high): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            [~,p,~,stat] = ttest2(results(3).x(c,2),results(3).x(~c,2));
            disp(['m_k (low vs. high): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            [~,p,~,stat] = ttest(results(3).x(:,1));
            disp(['m_alpha: t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            [~,p,~,stat] = ttest(results(3).x(:,2));
            disp(['m_k: t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            
            barerrorbar(m',se'); colormap summer;
            set(gca,'Fontsize',25,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'m_\alpha' 'm_k'});
            ylabel('Parameter estimate','FontSize',25);
            legend({'Low variance' 'High variance'},'FontSize',25,'Box','Off');
            
        case 'chavez_parameters'
            
            if nargin < 3; data = load_data('chavez'); end
            if nargin < 2; load results/results_chavez; end
            
            results(3).x(:,2) = -results(3).x(:,2); % flip sign on discounting magnitude effect
            
            m = mean(results(3).x(:,1:2));
            se = std(results(3).x(:,1:2))./sqrt(size(results(3).x,1));
            
            subplot(1,2,1);
            barerrorbar(m',se'); colormap summer
            set(gca,'Fontsize',25,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'m_\alpha' 'm_k'});
            ylabel('Parameter estimate','FontSize',25);
            [~,p,~,stat] = ttest(results(3).x(:,1:2));
            disp(['t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p)]);
            
            N = 7;
            [m,se] = quantile_stats(results(3).x(:,1),results(3).x(:,2),N);
            subplot(1,2,2);
            p = linspace(0,1,N);
            p = p(1:end-1) + diff(p)/2;
            errorbar(p,m,se,'-ok','LineWidth',4,'MarkerFaceColor','k','MarkerSize',12);
            xlabel('m_k quantile','FontSize',25);
            ylabel('m_\alpha','FontSize',25);
            set(gca,'FontSize',25,'XLim',[0 1]);
            
            set(gcf,'Position',[200 200 1200 450])
            
        case 'hardisty'
            
            r = [10 100 1000];
            beta = 0.2; s2_u = 10;
            for i = 1:2
                if i==1
                    s2_e = s2_u./(r.*beta);
                else
                    s2_e = s2_u./(-r.*beta);
                end
                k(i,:) = s2_e./s2_u;
            end
            bar(k);
            set(gca,'YLim',[min(k(:))-std(k(:))./2 max(k(:))+std(k(:))/2],'FontSize',25,'XTickLabel',{'Gain' 'Loss'});
            xlabel('Reward magnitude','FontSize',25);
            ylabel('Discount factor (k)','FontSize',25);
            legend({'$10' '$100' '$1000'},'FontSize',25,'Box','Off');
            
        case 'ballard'
            
            beta_inv = linspace(1e-10,1,4);
            for i = 1:length(beta_inv)
                s2_u = 1;
                r = [20+0.25*20 20; 2000+0.025*2000 2000];
                D = [14 0; 14 0];  % delays
                s2_e = beta_inv(i)*s2_u./r;
                k = s2_e./s2_u;
                d = 1./(1 + (k.*D));
                alpha(i,:) = 1./sqrt(sum(((d.^2).*s2_e.*D),2));
            end
            
            clr = [0 0 0; 0.5 0.5 0.5];
            for c = 1:2
                plot(beta_inv,log(alpha(:,c)),'-o','LineWidth',4,'MarkerSize',10,'Color',clr(c,:),'MarkerFaceColor',clr(c,:));
                hold on;
            end
            set(gca,'FontSize',25,'XLim',[min(beta_inv)-std(beta_inv)/2 max(beta_inv)+std(beta_inv)/2],'YLim',[0 20]);
            legend({'Low reward' 'High reward'},'FontSize',25,'box','off');
            xlabel('Attentional cost (\beta^{-1})','FontSize',25);
            ylabel('Inverse temperature (log \alpha)','FontSize',25);
            
            
        case 'schematic'
            
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
            
    end
    
end

function [m,se,X] = quantile_stats(x,y,N)
    
    q = quantile(y,N);
    
    for i = 1:length(q)-1
        ix = y>q(i) & y<=q(i+1);
        X{i} = x(ix);
        m(i) = nanmean(x(ix));
        se(i) = nanstd(x(ix))./sqrt(sum(~isnan(x(ix))));
    end
end

function varargout=shadedErrorBar(x,y,errBar,varargin)
    % generate continuous error bar area around a line plot
    %
    % function H=shadedErrorBar(x,y,errBar, ...)
    %
    % Purpose
    % Makes a 2-d line plot with a pretty shaded error bar made
    % using patch. Error bar color is chosen automatically.
    %
    %
    % Inputs (required)
    % x - vector of x values [optional, can be left empty]
    % y - vector of y values or a matrix of n observations by m cases
    %     where m has length(x);
    % errBar - if a vector we draw symmetric errorbars. If it has a size
    %          of [2,length(x)] then we draw asymmetric error bars with
    %          row 1 being the upper bar and row 2 being the lower bar
    %          (with respect to y -- see demo). ** alternatively **
    %          errBar can be a cellArray of two function handles. The
    %          first defines statistic the line should be and the second
    %          defines the error bar.
    %
    % Inputs (optional, param/value pairs)
    % 'lineProps' - ['-k' by default] defines the properties of
    %             the data line. e.g.:
    %             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
    % 'transparent' - [true  by default] if true, the shaded error
    %               bar is made transparent. However, for a transparent
    %               vector image you will need to save as PDF, not EPS,
    %               and set the figure renderer to "painters". An EPS
    %               will only be transparent if you set the renderer
    %               to OpenGL, however this makes a raster image.
    % 'patchSaturation'- [0.2 by default] The saturation of the patch color.
    %
    %
    %
    % Outputs
    % H - a structure of handles to the generated plot objects.
    %
    %
    % Examples:
    % y=randn(30,80);
    % x=1:size(y,2);
    %
    % 1)
    % shadedErrorBar(x,mean(y,1),std(y),'lineprops','g');
    %
    % 2)
    % shadedErrorBar(x,y,{@median,@std},'lineprops',{'r-o','markerfacecolor','r'});
    %
    % 3)
    % shadedErrorBar([],y,{@median,@(x) std(x)*1.96},'lineprops',{'r-o','markerfacecolor','k'});
    %
    % 4)
    % Overlay two transparent lines:
    % clf
    % y=randn(30,80)*10;
    % x=(1:size(y,2))-40;
    % shadedErrorBar(x,y,{@mean,@std},'lineprops','-r','transparent',1);
    % hold on
    % y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
    % shadedErrorBar(x,y,{@mean,@std},'lineprops','-b','transparent',1);
    % hold off
    %
    %
    % Rob Campbell - November 2009
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parse input arguments
    narginchk(3,inf)
    
    params = inputParser;
    params.CaseSensitive = false;
    params.addParameter('lineProps', '-k', @(x) ischar(x) | iscell(x));
    params.addParameter('transparent', true, @(x) islogical(x) || x==0 || x==1);
    params.addParameter('patchSaturation', 0.2, @(x) isnumeric(x) && x>=0 && x<=1);
    
    params.parse(varargin{:});
    
    %Extract values from the inputParser
    lineProps =  params.Results.lineProps;
    transparent =  params.Results.transparent;
    patchSaturation = params.Results.patchSaturation;
    
    if ~iscell(lineProps), lineProps={lineProps}; end
    
    
    %Process y using function handles if needed to make the error bar dynamically
    if iscell(errBar)
        fun1=errBar{1};
        fun2=errBar{2};
        errBar=fun2(y);
        y=fun1(y);
    else
        y=y(:).';
    end
    
    if isempty(x)
        x=1:length(y);
    else
        x=x(:).';
    end
    
    
    %Make upper and lower error bars if only one was specified
    if length(errBar)==length(errBar(:))
        errBar=repmat(errBar(:)',2,1);
    else
        s=size(errBar);
        f=find(s==2);
        if isempty(f), error('errBar has the wrong size'), end
        if f==2, errBar=errBar'; end
    end
    
    if length(x) ~= length(errBar)
        error('length(x) must equal length(errBar)')
    end
    
    
    %Log the hold status so we don't change
    initialHoldStatus=ishold;
    if ~initialHoldStatus, hold on,  end
    
    H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation);
    
    if ~initialHoldStatus, hold off, end
    
    if nargout==1
        varargout{1}=H;
    end
    
end

function H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot to get the parameters of the line
    
    H.mainLine=plot(x,y,lineProps{:},'LineWidth',4);
    
    
    % Work out the color of the shaded region and associated lines.
    % Here we have the option of choosing alpha or a de-saturated
    % solid colour for the patch surface.
    mainLineColor=get(H.mainLine,'color');
    edgeColor=mainLineColor+(1-mainLineColor)*0.55;
    
    if transparent
        faceAlpha=patchSaturation;
        patchColor=mainLineColor;
    else
        faceAlpha=1;
        patchColor=mainLineColor+(1-mainLineColor)*(1-patchSaturation);
    end
    
    
    %Calculate the error bars
    uE=y+errBar(1,:);
    lE=y-errBar(2,:);
    
    
    %Add the patch error bar
    
    
    
    %Make the patch
    yP=[lE,fliplr(uE)];
    xP=[x,fliplr(x)];
    
    %remove nans otherwise patch won't work
    xP(isnan(yP))=[];
    yP(isnan(yP))=[];
    
    
    if(isdatetime(x))
        H.patch=patch(datenum(xP),yP,1,'HandleVisibility','off');
    else
        H.patch=patch(xP,yP,1,'HandleVisibility','off');
    end
    
    set(H.patch,'facecolor',patchColor, ...
        'edgecolor','none', ...
        'facealpha',faceAlpha)
    
    
    %Make pretty edges around the patch.
    H.edge(1)=plot(x,lE,'-','color',edgeColor,'HandleVisibility','off');
    H.edge(2)=plot(x,uE,'-','color',edgeColor,'HandleVisibility','off');
    
    
    
    uistack(H.mainLine,'top') % Bring the main line to the top
end