function [my_fig,plot_params] = multi_run_GUI(filespath,xaxis,tile,ln_clr,aggregate,fixed_params)

persistent prev_filespath params vals my_efficient_map_tmp dims_cumprod;
if nargin==0
    params.workers=0;
    [~,~,~,~,params,filespath] = multi_run_calculate(params,true);
    xaxis = "graph_tree";
    ln_clr = "embedding_method";
    tile = "clustered_entity";
end

if ~exist('filespath','var')
    error("Not enough inputs.")
end

% TODO: Make these initializations even when variables exist but are empty
if ~exist('xaxis','var')
    xaxis = "";
end
if ~exist('tile','var')
    tile = [""];
end
if ~exist('ln_clr','var')
    ln_clr = [""];
end
if ~exist('aggregate','var')
    % TODO: allow this to be empty
    aggregate = "sample";
end

if ~exist('fixed_params','var')
    fixed_params = [];
end

if isempty(vals) || isempty(my_efficient_map_tmp) || isempty(dims_cumprod) || isempty(params) || ~isequal(filespath,prev_filespath)
    tic
    prev_filespath = filespath;
    filenames = string({dir(filespath).name})';
    cluster_filenames = filenames(contains(filenames,"clusterings_part"));
    clusterings = [];
    for i = 1:numel(cluster_filenames)
        % TODO: Introduce a loading bar
        disp("Loading data: "+round(i/numel(cluster_filenames)*100)+"%")
        l = load(filespath+"/"+cluster_filenames(i));
        if isfield(l,'params')
            params = l.params;
        else
            clusterings = [clusterings l.clusterings];
        end
    end
    toc
    if any(filenames=="efficient_map_data.mat")
        disp("Loading Efficient Map Data: Started")
        tic
        l = load(filespath+"/efficient_map_data.mat");
        linear_inds = l.linear_inds;
        dims_cumprod = l.dims_cumprod;
        my_efficient_map_size = l.my_efficient_map_size;
        toc
        disp("Loading Efficient Map Data: Finished")
        disp("Creating Efficient Map: Started")
        tic
    else
        disp("Creating Efficient Map: Started")
        tic
        dims_num = max(cellfun(@numel,{clusterings.params_inds}));
        tmp = cellfun(@(x) [x , zeros(1,dims_num-numel(x))],{clusterings.params_inds},'UniformOutput',false);
        params_inds_all = cat(1,tmp{:})+1;
        dim_sizes = max(params_inds_all,[],1);
        dims_cumprod = cumprod([1 dim_sizes(1:end-1)]);
        linear_inds_tmp = (params_inds_all-1)*dims_cumprod'+1;
        [linear_inds(:,1),linear_inds(:,2)] = ind2sub([2^48-1, 2^48-1],linear_inds_tmp);
        my_efficient_map_size = [min(2^48-1,max(linear_inds_tmp)),max(linear_inds(:,2))] ;
    end
    my_efficient_map_tmp = sparse(linear_inds(:,1),linear_inds(:,2), ...
        [1:size(linear_inds,1)]',my_efficient_map_size(1),my_efficient_map_size(2));
    vals = {clusterings.data};
    clear clusterings
    toc
    disp("Creating Efficient Map: Finished")
    clear l
    if ~any(filenames=="efficient_map_data.mat")
        disp("Saving Efficient Map: Started")
        tic
        save(filespath+"/efficient_map_data.mat",'linear_inds','dims_cumprod','my_efficient_map_size','-v7.3');
        toc
        disp("Saving Efficient Map: Finished")
    end
    clear linear_inds 
else
    disp("Using loaded data")
end
disp("~~~~~~~")

    function value = my_efficient_map(key)
        tmp = key*dims_cumprod(1:numel(key))'+1;
        % TODO: 2^48-1 is for 64 bit systems only and could change in future matlab versions
        if tmp <= 2^48-1
            value = vals{my_efficient_map_tmp(tmp)};
        else
            %             [tmp1,tmp2] = ind2sub(size(my_efficient_map_tmp),tmp);
            tmp1 = mod(tmp-1,2^48-1)+1;
            tmp2 = ceil(tmp/(2^48-1));
            value = vals{my_efficient_map_tmp(tmp1,tmp2)};
        end
    end



prms = params;
for i = 1:numel(fixed_params)
    evalc("prms."+fixed_params(i));
end

% Keeps only the first clustering measure and entity to avoid comparisons
% between different types of them. If any of these two has been assigned to
% some dimension of the visualization, they are set to the desired value as
% needed
if numel(prms.clustering_measure)>1
    prms.clustering_measure = params.clustering_measure(1);
end
if numel(prms.clustered_entity)>1
    prms.clustered_entity = params.clustered_entity(1);
end

tic_all = tic;
% TODO: check whether aggregates need to be added to the list

% If a parameter which is not at the first level of the parameter tree is
% selected to be scanned, then all of its parent parameters in the previous
% tree levels are fixed so that no sibling of any parent is considered in
% any of the calculations.
dims_names = ["ln_clr","tile","xaxis"];
clustering_measure_is_scanned = false;
clustered_entity_is_scanned = false;
cur_dim = [];
for i = 1:numel(dims_names)
    evalc("cur_dim = "+dims_names(i));
    if ~isequal(cur_dim,"")
        for k = 1:numel(cur_dim)
            if cur_dim(k)=="clustering_measure"
                clustering_measure_is_scanned = true;
            elseif cur_dim(k)=="clustered_entity"
                clustered_entity_is_scanned = true;
            end
            cur_dim_parts = strsplit(cur_dim(k),'.');
            tmp = "prms";
            for j = 1:2:(numel(cur_dim_parts)-1)
                tmp = tmp+"."+cur_dim_parts{j};
                evalc("tmp2="+tmp);
                if isstruct(tmp2)
                    evalc(tmp+".choice = string('"+cur_dim_parts{j+1}+"')");
                end
                tmp = tmp+"."+cur_dim_parts{j+1};
            end
        end
    end
end



xaxis_dat =  eval("params."+xaxis);
if isstruct(xaxis_dat)
    xaxis_all = string(fieldnames(xaxis_dat));
else
    xaxis_all = xaxis_dat ;
end

% Combine the following two if statements 
if ~isequal(ln_clr,"")
    [ln_clr_all,ln_clr_dat] = deal(cell(0));
    ln_clr_sizes = [];
    for i = 1:numel(ln_clr)
        ln_clr_dat{i} =  eval("params."+ln_clr(i));
        if isstruct(ln_clr_dat{i})
            ln_clr_all{i} = string(fieldnames(ln_clr_dat{i}));
        else
            ln_clr_all{i} = ln_clr_dat{i};
        end
        ln_clr_sizes(i) = numel(ln_clr_all{i});
    end
else
    ln_clr_sizes = 1;
end

if ~isequal(tile,"")
    [tile_all,tile_dat] = deal(cell(0));
    tile_sizes = [];
    for i = 1:numel(tile)
        tile_dat{i} =  eval("params."+tile(i));
        if isstruct(tile_dat{i})
            tile_all{i} = string(fieldnames(tile_dat{i}));
        else
            tile_all{i} = tile_dat{i};
        end
        tile_sizes(i) = numel(tile_all{i});
    end
else
    tile_sizes = 1;
end


% TODO: Aggregate over multiple variables with multiple types of aggregation functions
if numel(aggregate)>1
    error("Only one aggregate is allowed")
end
aggregate_dat =  eval("prms."+aggregate);
if isstruct(aggregate_dat)
    error("The children of the aggregate must be leaves in the parameter tree.");
    aggregate_all = string(fieldnames(aggregate_dat));
else
    aggregate_all = aggregate_dat;
end


plot_params = cell(prod(tile_sizes),prod(ln_clr_sizes),numel(xaxis_all));
my_colors =  [
    0       0.4470  0.7410
    0.8500  0.3250  0.0980
    0.9290  0.6940  0.1250
    0.4940  0.1840  0.5560
    0.4660  0.6740  0.1880
    0.3010  0.7450  0.9330
    0.6350  0.0780  0.1840
    0       0       1
    0       0.5     0
    1       0       0
    0       0.75    0.75
    0.75    0       0.75
    0.75    0.75    0
    0.25    0.25    0.25];
my_markers = {'+','o','*','.','x','square','diamond','v','^','>','<','pentagram','hexagram'};
% close all
my_fig = figure;
if numel(tile_sizes) ~= 2
    tiledlayout(my_fig,'flow','TileSpacing', 'none', 'Padding', 'none')
else
    tiledlayout(my_fig,tile_sizes(2),tile_sizes(1),'TileSpacing', 'none', 'Padding', 'none')
end


points_plotted = 0;
total_points_to_plot = prod([prod(tile_sizes), prod(ln_clr_sizes), numel(xaxis_all)]);
for tile_ind = 1:prod(tile_sizes)
    nexttile(tile_ind)
    if ~clustering_measure_is_scanned
        ylabel(strrep(prms.clustering_measure,'_','\_'));
    end  
    cur_title = strings(0);
    if ~isequal(tile,"")
        tile_subind = cell(1,numel(tile));
        [tile_subind{:}] = ind2sub(tile_sizes,tile_ind);
        cur_title = strings(1,numel(tile));
        j = 0;
        for i = 1:numel(tile)
            if isstruct(tile_dat{i})
                evalc("prms."+tile(i)+".choice=tile_all{i}(tile_subind{i})");
            else
                evalc("prms."+tile(i)+"=tile_all{i}(tile_subind{i})");
            end
            try
                cur_title(i) = tile(i) + ": " +string(tile_all{i}(tile_subind{i}));
            catch
                cur_title(i) = tile(i) + ": " +  string(tile_subind{i});
            end
            cur_title(i) =  strrep(cur_title(i),'_','\_');
            if tile(i)=="clustering_measure" 
                ylabel(strrep(tile_all{i}(tile_subind{i}),'_','\_'));
                j=i;
            end
        end
        title(cur_title([1:j-1,j+1:end]))
    end
    hold on
    cur_legend = strings(1,prod(ln_clr_sizes));
    for ln_clr_ind = 1:prod(ln_clr_sizes)
        if ~isequal(ln_clr,"")
            ln_clr_subind = cell(1,numel(ln_clr));
            [ln_clr_subind{:}] = ind2sub(ln_clr_sizes,ln_clr_ind);
            for i = 1:numel(ln_clr)
                if isstruct(ln_clr_dat{i})
                    evalc("prms."+ln_clr(i)+".choice=ln_clr_all{i}(ln_clr_subind{i})");
                else
                    evalc("prms."+ln_clr(i)+"=ln_clr_all{i}(ln_clr_subind{i})");
                end

                try
                    tmp = strrep( string(ln_clr_all{i}(ln_clr_subind{i})),'_','\_');
                catch
                    tmp = string(ln_clr_subind{i});
                end
                cur_legend(ln_clr_ind) = cur_legend(ln_clr_ind)+tmp;
                if i < numel(ln_clr)
                    cur_legend(ln_clr_ind) = cur_legend(ln_clr_ind)+", ";
                end
            end
        end
        data = [];
        
        ft = [];
        
        for xaxis_ind = 1:numel(xaxis_all)
            if isstruct(xaxis_dat)
                evalc("prms."+xaxis+".choice=xaxis_all(xaxis_ind)");
            else
                evalc("prms."+xaxis+"=xaxis_all(xaxis_ind)");
            end
            
            tic_cur = tic;

            % TODO: selected a random sample to generate all combinations 
            evalc("prms."+aggregate+" = aggregate_all(1)");
            prms2 = prune_values(prms);
            prms2_all = generate_combinations(prms2);
            data_eff_inds = cell(numel(aggregate_all),numel(prms2_all));
            parfor(aggregate_2_ind = 1:numel(prms2_all),params.workers)
                prms2_all_parfor = prms2_all;
                data2_eff_inds = cell(1,numel(aggregate_all)); 
                data_tmp = [];
                [cur_inds,cur_names] = params2inds(prms2_all_parfor{aggregate_2_ind}.params,params);

                % TODO: make it more flexible
                aggregate_ind_ind = find(cur_names==aggregate);
                
                for aggregate_ind = 1:numel(aggregate_all)
                    cur_inds(aggregate_ind_ind) = aggregate_ind;
                    data2_eff_inds{aggregate_ind} = cur_inds;
                end
                data_eff_inds(:,aggregate_2_ind) = data2_eff_inds;
                %                 data(:,aggregate_2_ind) = data_tmp;
            end

            data = cellfun(@(x)mean(my_efficient_map(x)),data_eff_inds);
            ft_tmp = prctile(data,[25 50 75],1);
            % TODO: This assumes that the maximum value corresponds to the best performance. Make it more flexible.
            % TODO: Compares the performances of various parameter combinations based on their median performance. Make it more flexible.
            [~,ft_best_ind] = max(ft_tmp(2,:),[],2);
            ft(:,xaxis_ind) = ft_tmp(:,ft_best_ind);
            
            % TODO: Check soundness
            [~,best_aggregate_ind] = min(abs(ft_tmp(2,ft_best_ind)-data(:,ft_best_ind)));
            prms2_all{ft_best_ind}.params.sample = aggregate_all(best_aggregate_ind);
            plot_params(tile_ind,ln_clr_ind,xaxis_ind) = prms2_all(ft_best_ind);
            
            points_plotted = points_plotted+1;
            
            tmp = round(points_plotted/total_points_to_plot*100);
            tmp2 = "Creating figure: "+tmp+...
                "% | Iteration time: "+char(seconds(toc(tic_cur)),'hh:mm:ss')+...
                " | Total time: "+char(seconds(toc(tic_all)),'hh:mm:ss')+...
                " | TILE: "+strjoin(cur_title,", ")+...
                " | LINE: "+strjoin(ln_clr,", ")+": "+cur_legend(ln_clr_ind);
            try
                tmp2 = tmp2+ " | XAXIS: "+xaxis+": "+xaxis_all(xaxis_ind);
            catch
                tmp2 = tmp2+ " | XAXIS: "+xaxis+": "+xaxis_ind;
            end
            
            disp(strrep(tmp2,"\_","_"));
            
        end
        
        if isnumeric(xaxis_all)
            x_vals = xaxis_all;
        else
            x_vals = 1:numel(xaxis_all);
        end
        
        
        if numel(ln_clr_sizes) ~= 2
            fill([x_vals';flipud(x_vals')],[ft(3,:)';flipud(ft(1,:)')],my_colors(ln_clr_ind,:),'FaceAlpha',0.15,'linestyle','none','HandleVisibility','off');
            plot(x_vals,ft(2,:),'Color',my_colors(ln_clr_ind,:),'Linewidth',1);
        else
            fill([x_vals';flipud(x_vals')],[ft(3,:)';flipud(ft(1,:)')],my_colors(ln_clr_subind{1},:),'FaceAlpha',0.15,'linestyle','none','HandleVisibility','off');
            plot(x_vals,ft(2,:),'Marker',my_markers{ln_clr_subind{2}},'Color',my_colors(ln_clr_subind{1},:),'Linewidth',1);
        end
        xlabel(strrep(xaxis,'_','\_'))
        drawnow
    end
    if ~isnumeric(xaxis_all)
        xticks(nexttile(tile_ind),1:numel(xaxis_all));
        try
            xticklabels(nexttile(tile_ind),string(xaxis_all')');
        catch
        end
    end
    hold off
    grid on
    if ~isequal(ln_clr,"")
        lgd = legend(cur_legend,'Location','Best');
        lgd.Title.String = strrep(ln_clr,'_','\_');
        lgd.Title.FontSize = 6;
        lgd.EdgeColor  = "none";
        set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.2;.2;.2;.0])) ;
%         lgd.Box = 'off';
    end
    axis tight
end
linkaxes(my_fig.Children.Children)

%set(findobj('Type','legend'),'Box','off')

disp("Total time:" + char(seconds(toc(tic_all)),'hh:mm:ss'))
end

function prms = prune_values(prms)
% Removes any parameter values that cannot coexist with the values of the
% parameters for which a value has been selected.
f = fieldnames(prms);
for i = 1:numel(f)
    if isstruct(prms.(f{i}))
        f2 = fieldnames(prms.(f{i}));
        if isfield(prms.(f{i}),'choice')
            cur_choice = prms.(f{i}).choice;
        else
            cur_choice = [];
        end

        for j = 1:numel(f2)
            if ~isempty(cur_choice) && ~isequal(f2{j},cur_choice)
                prms.(f{i}) = rmfield(prms.(f{i}),f2{j});
            elseif isstruct(prms.(f{i}).(f2{j}))
                prms.(f{i}).(f2{j}) = prune_values(prms.(f{i}).(f2{j}));
            end
        end
    end
end

end



% ~~~~~~~~~
% set(findobj('Type','legend'),'Visible','on')
% ~~~~~~~~~
%
% aggFun_list = {@(x)median(x,1),@(x)mean(x,1),@(x)max(x,[],1),@(x)min(x,[],1)};
% figure
% for figures_ind=1:size(data.U_all,3)
%     nexttile(figures_ind)
%     for lines_ind=1:size(data.U_all,2)
%         nmi_avg_all=[];
%         obj_all=[];
%         for x_ticks_ind=1:size(data.U_all,1)
%            parfor aggregates_ind=1:size(data.U_all,4)
%                 U = data.U_all{x_ticks_ind,lines_ind,figures_ind,aggregates_ind};
%                 A = data.A_all{x_ticks_ind,lines_ind,figures_ind,aggregates_ind};
%                 B = data.B_all{x_ticks_ind,lines_ind,figures_ind,aggregates_ind};
%                 nodes_labels = data.nodes_labels_all{x_ticks_ind,lines_ind,figures_ind,aggregates_ind};
%                 views_labels = data.views_labels_all{x_ticks_ind,lines_ind,figures_ind,aggregates_ind};
%
%
%                 U_nrm_type_ind=[3 2];
%                 U_clus_type_ind=1;
%                 %                 U_cl_num = -ones(1,numel(params.graph_tree.Children));
%                 U_cl_num = [3 2 2];
%
%                 pred = cluster_embeddings(U,A,B,nodes_labels,views_labels,U_nrm_type_ind,U_clus_type_ind,U_cl_num);
%
%                 U_cl_qual_all(aggregates_ind,x_ticks_ind) = pred.nodes.cluster_qual.silhouette;
%                 obj_all(aggregates_ind,x_ticks_ind) = data.obj_cur_all{x_ticks_ind,lines_ind,figures_ind,aggregates_ind};
%                 nmi_avg_all(aggregates_ind,x_ticks_ind) = pred.nodes.cluster_qual.NMI;
%             end
%         end
% %         ~~~~~~~~~~ finds nmi_avg with smallest objective function value
% %                 [~,min_obj_inds]= min(obj_all,[],1);
% %                 inds = sub2ind(size(nmi_avg_all),min_obj_inds,[1:size(nmi_avg_all,2)]);
% %                 nmi_avg_all = nmi_avg_all(inds);
%
% %         ~~~~~~~~~~ finds nmi_avg with largest silhouette value
% %                 [~,max_U_cl_qual_inds]= max(U_cl_qual_all,[],1);
% %                 inds = sub2ind(size(nmi_avg_all),max_U_cl_qual_inds,[1:size(nmi_avg_all,2)]);
% %                 nmi_avg_all = nmi_avg_all(inds);
%
%         % ~~~~~~~~~~ calculates median mean max min
% %                 nmi_avg_all = aggFun_list{3}(nmi_avg_all);
%
% %                 plot(x_vals,nmi_avg_all)
% %                 hold on
%
%         ft = prctile(nmi_avg_all,[25 50 75],1);
%         fill([x_vals';flipud(x_vals')],[ft(3,:)';flipud(ft(1,:)')],my_colors(lines_ind,:),'FaceAlpha',0.15,'linestyle','none','HandleVisibility','off');
%         hold on
%         plot(x_vals,ft(2,:),'Color',my_colors(lines_ind,:),'Linewidth',1);
%
%         drawnow
%
%     end
%     hold off
%
%     ylim([0 1])
%     grid on
%
%     if ~isempty(x_ticks)
%         xticks(x_vals)
%         xticklabels(replace(x_ticks,"_","\_"))
%     end
%
%     if ~isempty(figures)
%         title(replace(figures(figures_ind),"_","\_"))
%     end
%     xlabel(replace(x_label,"_","\_"))
%     if ~isempty(lines)
%         legend(replace(lines,"_","\_"),'Location','Best')
%     end
% end
