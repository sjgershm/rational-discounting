function h = mytitle(txt,just,varargin)
    
    % Plot titles with left or right justification.
    %
    % USAGE: h = mytitle(txt,[just],[varargin])
    %
    % INPUTS:
    %   txt - title string
    %   just (optional) - justification: 'Left' (default) or 'Right'
    %   varargin (optional) - additional arguments for title
    %
    % OUTPUTS:
    %   h - handle to title object
    %
    % Sam Gershman, Dec 2011
    
    if nargin < 2 || isempty(just)
        just = 'Left';
    end
    
    switch just
        case 'Left'
            h = title(txt,'HorizontalAlignment','left','Units','normalized','Position',[0 1],varargin{:});
        case 'Right'
            h = title(txt,'HorizontalAlignment','right','Units','normalized','Position',[1 1],varargin{:});
    end