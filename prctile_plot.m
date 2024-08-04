function [line_handle,patch_handle] = prctile_plot(varargin)
    %PRCTILE_PLOT Calculates and plots a series of simplified boxplots.
    %
    %    PRCTILE_PLOT(x,Y) calculates and plots the median and the 25-th 
    %    and 75-th percentiles of each column of Y versus the corresponding
    %    value of x.
    %
    %    PRCTILE_PLOT(Y) plots the same quantities but versus the indices 
    %    of Y.
    %    
    %    PRCTILE_PLOT(ax,...) plots into the axes with handle ax.
    %
    %    PRCTILE_PLOT(...,Name,Value) specifies plot properties using one
    %    or more Name,Value pair arguments. Any pair compatible with the 
    %    built-in PLOT function can be used to customize the line 
    %    representing the medians. Additionally, the 'FaceAlpha' property 
    %    can be used to customize the transparency of the area representing
    %    the 25-th and 75-th percentiles.  
    %    
    %    [line_handle,patch_handle] = PRCTILE_PLOT(...) returns handles to
    %    the plotted line and patch objects.

    if nargin == 0
        error("Not enough input arguments.")
    else
        rest_args = [];
        if isgraphics(varargin{1},'axes')
            ax = varargin{1};
            rest_args = varargin(2:end);
        else
            rest_args = varargin;
        end
        if numel(rest_args)==0
            error("Not enough input arguments.")
        else
            if ~isnumeric(rest_args{1})
                error("Input data must be a numeric array");
            elseif numel(rest_args)==1 || ~isnumeric(rest_args{2})
                Y = rest_args{1};
                rest_args = rest_args(2:end);
            else
                x = rest_args{1};
                Y = rest_args{2};
                rest_args = rest_args(3:end);
            end
            for prop = ["Color","FaceAlpha"]
                tmp = find(cellfun(@(x)isequal(x,prop),rest_args));
                if ~isempty(tmp)
                    if  numel(rest_args)==tmp 
                        error("Not enough input arguments.")
                    else
                        switch prop
                            case "Color"
                                color = rest_args{tmp+1};
                            case "FaceAlpha"
                                face_alpha = rest_args{tmp+1};
                        end
                        rest_args(tmp:tmp+1) = [];
                    end
                end
            end
        end
    end
    for cur_var_name = ["x","ax","color","face_alpha"] 
        if ~exist(cur_var_name,'var') || isempty(eval(cur_var_name))
            switch cur_var_name
                case 'x'
                    x = 1:size(Y,2);
                case 'ax'
                    ax = gca;
                case 'color'
                    colors_all = get(ax,'colororder');
                    if ishold(ax)
                        color = colors_all(mod(numel(ax.Children),size(colors_all,1))+1,:);
                    else
                        color = colors_all(1,:);
                    end
                case 'face_alpha'
                    face_alpha = 0.15;
            end
        end
    end
    hold_was_on = ishold(ax);
    try
        ft = prctile(Y,[25 50 75],1);
        patch_handle = fill(ax,[x';flipud(x')],[ft(3,:)';flipud(ft(1,:)')],color,'FaceAlpha',face_alpha,'linestyle','none','HandleVisibility','off');
        hold(ax,'on');
        line_handle = plot(ax,x,ft(2,:),'Color',color);
        while ~isempty(rest_args)
            if numel(rest_args)==1
                error("Not enough input arguments.")
            end
            set(line_handle,rest_args{1},rest_args{2});
            rest_args(1:2) = [];
        end
    catch 
        if ~hold_was_on
            hold(ax,'off');
        end
        rethrow(lasterror);
    end
    if ~hold_was_on
        hold(ax,'off');
    end
end
