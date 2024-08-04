function single_run_GUI(varargin)
%SINGLE_RUN_GUI gives detailed visualizations of a single multi-view graph
%clustering
%
%   SINGLE_RUN_GUI(graph_tree) takes as input a graph_tree_root object
%   which defines an artificial graph. The type and the size of each node
%   cluster along with the size of each view cluster should be defined.
%   For example, below we define 2 view clusters with 3 views each such that
%   the first view cluster has 2 cliques with sizes 20 and 10,
%   respectively, and the second view cluster has 2 cliques with sizes 5
%   and 25, respectively:
%       par_all = graph_tree_root;
%       par_all.Children(1).slices_num = 3;
%       par_all.Children(1).Children(1).type = 'clique';
%       par_all.Children(1).Children(1).size = 20; 
%       par_all.Children(1).Children(2).type = 'clique';
%       par_all.Children(1).Children(2).size = 10; 
%       par_all.Children(2).slices_num = 3;
%       par_all.Children(2).Children(1).type = 'clique';
%       par_all.Children(2).Children(1).size = 5; 
%       par_all.Children(2).Children(2).type = 'clique';
%       par_all.Children(2).Children(2).size = 25; 
%
%   SINGLE_RUN_GUI(X) takes as input the adjaceny tensor X such that
%   X(:,:,i) is the adjacency matrix of the i-th view.
%
%   SINGLE_RUN_GUI(X,view_labels,node_labels) considers view labels defined
%   in the row vector view_labels, and node labels of the i-th view cluster
%   in node_labels{i} as a row vector. Note that the i-th view cluster is
%   considered to be the one corresponding to the i-th value of the sorted
%   unique view labels as defined by unique(view_labels).
%   
%   Finally, there are the following commands:
%   .-Go to the next view
%   ,-Go to the previous view
%   '-Restrict navigation only to views corresponding the current view
%   cluster
%   /-Go to the next view cluster
%   Tab-Switch navigation type between ground-truth and calculated view
%   clusters
%   e-Select the next parameter
%   q-Select the previous parameter 
%   w-Select the next value for the parameter (increases value for
%   numerical parameters)
%   s-Select the previous value for the parameter (decreases value for
%   numerical parameters)
%   d-Select the next sub-parameter
%   a-Select the previous sub-parameter
%   Space-Recalculate everything based on the selected parameters
%   Enter-Recalculate only clustering (can also be used for refreshing plots)
%   h-Hide legends (removing legends improves responsiveness of GUI)
%   Press g to toggle the graph modification mode:
%       Enter-Select the first node cluster
%       Right arrow-Go to the next node cluster
%       Left arrow-Go to the previous node cluster
%       Up arrow-Assign the next community type to the current node cluster
%       Down arrow-Assign the previous community type to the current node cluster
%       +-Increase the size of the selected node cluster
%       --Decrease the size of the selected node cluster
%       i-Inseart new node cluster
%       d-Deleted selected node cluster
%       press n to toggle the noise and sparsity modification modes:
%           +-Increase the amount of noise/sparsity
%           --Decrease the amount of noise/sparsity
%       press
%       Press s to toggle the view clustering modification mode:
%           +-Insert new view cluster
%           --Remove view cluster corresponding to the currently selected view
%           i-Inseart new view cluster
%           d-Deleted selected view cluster

% TODO: Should be able to read precalculated data from file.
% TODO: Order clusters in Clustered Adjacency Matrix as in Original Adjacency matrix
% FIXME: NMI fails after data are modified by removing a cluster
% FIXME: predicted view labels in legends do not match the exact ID of the
% corresponding ground-truth view label.
    % -----------------------------------------
    % ------------ Initializations ------------
    % -----------------------------------------
    mtimesx_exists = exist('mtimesx')==3;
    X = [];
    nodes_labels = [];
    nodes_labels_names = [];
    views_labels_names = [];
    views_labels = [];
    gui_data = [];
    init_par = [];
    cur_par = [];
    % rng(700)
    L_type_ind = 2;
    alg = 2;
    thres = 1e-6;
    max_iters = 1000;
    print_type = "basic";
    alg_opts = cell(1,3);
    alg_opts{2} = {0.01,0.00,1e-6}; % ComClus: beta, rho, thres_inner
    alg_opts{3} = {0,'random'}; % Symmetric Richcom: rho
    alg_opts{4} = {1,'random'}; % CMNC: delta
    U_nrm_type_ind  = [2 1]; % [1:skip  2:U*sqrt*(abs(diag(B)))  3: U*diag(B), 1:skip  2:unit vec]
    U_clus_type_ind = 3; % 1:k-means, 2: max abs value, 3: large inner products
    cl_num_auto = true;
    U_cl_num = [];

    k_orig = []; 
    clustered_matrix_handle = [];
    if nargin==1
        if isa(varargin{1},'graph_tree_root')
            par_all = varargin{1};
        elseif isa(varargin{1},'double')
            par_all = graph_tree_root;
            par_all.Data = varargin{1};
        else
            error('Invalid input.')
        end
    elseif nargin==3
        par_all = graph_tree_root;
        par_all.Data = varargin{1};
        par_all.labels = varargin{2};
        for i = 1:numel(varargin{3})
            par_all.Children(i).labels = varargin{3}{i};
        end
    elseif nargin~=0
        error('Invalid number of inputs')
    end
    graph_regen();

    %different initializations based on whether ground-truth labels exist
    if ~isempty(nodes_labels)
        for i=1:numel(nodes_labels)
            U_cl_num(i) = max([numel(unique(nodes_labels{i})), 1]);
        end
        if isempty(U_cl_num)
            U_cl_num = 1;
        end
        M = numel(U_cl_num);
    else
        %     com_num = 1;
        %     U_cl_num = 2*ones(1,com_num);
        M = 3;
        U_cl_num = [50 50 50];
    end
    R = sum(U_cl_num);
    largest_inner_prod_thres = 0.89;
    largest_inner_prod_step = 0.01; 

    while true
        if isfield(gui_data,'txt_handles')
            tmp = gui_data.status_txt_order_id;
            gui_data.txt_handles(tmp).String = "Calculating";
            drawnow
        end


        [U_nrm,eigvals]=deal([]);

        if ~isfield(gui_data,'force_plot_is_on') || ~gui_data.force_plot_is_on
            alg_opts_struct{2}.beta = alg_opts{2}{1};
            alg_opts_struct{2}.rho = alg_opts{2}{2};
            alg_opts_struct{2}.thres_inner = alg_opts{2}{3};
            alg_opts_struct{3}.rho = alg_opts{3}{1};
            alg_opts_struct{3}.structure = alg_opts{3}{2};
            alg_opts_struct{4}.delta = alg_opts{4}{1};
            alg_opts_struct{4}.structure = alg_opts{4}{2};
            [U,A,B,L,R,M,method_specific_vars,iters,fit_time,obj_cur] =...
                generate_embeddings(X,nodes_labels,alg,L_type_ind,R,M,...
                thres,max_iters,alg_opts_struct{alg},print_type);

            L_rec = parafac2full(U,U,A*B,mtimesx_exists);

            disp("fit: "+fit_time)
            if alg == 3
                U_cl_num(M+1:end)=[];
            end
        end

        % For estimating the number of communities in each view structure automatically
        if cl_num_auto
            U_cl_num(:) = -1;
        end

        tic

        % TODO: Update params code structure in entire GUI code.
        cur_params=[];
        switch U_clus_type_ind
            case 1
                cur_params.clustering_method.kmeans.replicates = 15;
                cur_params.clustering_method.kmeans.clusters_num = strjoin(string(U_cl_num)," ");
                switch U_nrm_type_ind(2)
                    case 1
                        cur_params.clustering_method.kmeans.row_normalization_type ="none";
                    case 2
                        cur_params.clustering_method.kmeans.row_normalization_type ="unit";
                end
            case 2
                cur_params.clustering_method.maximum.nofield = "nofield";
            case 3 % TODO: add in GUI
                cur_params.clustering_method.large_inner_prod.thres = largest_inner_prod_thres;
        end

        switch U_nrm_type_ind(1)
            case 1
                cur_params.column_normalization_type = "none";
            case 2
                cur_params.column_normalization_type = "sqrtB";
            case 3
                cur_params.column_normalization_type = "B";
        end
        cur_params.graph_tree = par_all;

        [pred,U_nrm] = cluster_embeddings(U,A,B,cur_params);
        disp("cluster: "+toc)

        if ~isfield(gui_data,'force_plot_is_on') || ~gui_data.force_plot_is_on
            selected_param_scan_name = "NMI";
            if ~isfield(pred.nodes.cluster_qual,selected_param_scan_name)
                selected_param_scan_name = 'silhouette_equal';
            end
            nodes_scan_qual = [];
            views_scan_qual = [];
            NMI_all = [];
            thres_all =  linspace(0.5,0.9999,10);
            thres_ind_all =  1:numel(thres_all);
            parfor sample = 1:0
                cur_params_2 = cur_params;
                [U_tmp,A_tmp,B_tmp] = generate_embeddings(X,nodes_labels,alg,L_type_ind,R,M,thres,max_iters,alg_opts_struct{alg},'nothing');
                for thres_ind = thres_ind_all
                    cur_params_2.clustering_method.large_inner_prod.thres = thres_all(thres_ind);
                    pred_2 = cluster_embeddings(U_tmp,A_tmp,B_tmp,cur_params_2);
                    nodes_scan_qual(sample,thres_ind,:) = ...
                        pred_2.nodes.cluster_qual.(selected_param_scan_name);
                    views_scan_qual(sample,thres_ind) = ...
                        pred_2.views.cluster_qual.(selected_param_scan_name);
                    %                     NMI_all(sample,thres_ind,:) = pred_2.nodes.cluster_qual.NMI;

                end 
            end
        end

        U_cl_num = pred.nodes.clusters_num;
        uniq_pred_view_labels = unique(pred.views.labels);
        perm_proj = cell(1,numel(uniq_pred_view_labels));
        for i = 1:numel(uniq_pred_view_labels)
            % Soft-matching node clusters to original labels
            j = uniq_pred_view_labels(i);
            if j <= numel(nodes_labels) && numel(unique(nodes_labels{j}))>1
                perm_proj{i} = (...
                    normalize_fibers(pred.nodes.labels{i}==unique(pred.nodes.labels{i})',1)'*...
                    normalize_fibers(nodes_labels{j}'==unique(nodes_labels{j}),1)...
                    ).^2;
                [~,tmp] = max(perm_proj{i},[],2);

                % TODO: let the user activate/deactivate the line below
                perm_proj{i} = tmp==[1:size(perm_proj{i},2)];
            else
                perm_proj{i} = eye(U_cl_num(i));
            end
        end


        %         A,B
        gui_plot()
        gui_data.force_plot_is_on = false;
        gui_data.recalculate  = false;
        if ~gui_data.navigate_based_on_orig_structure
            gui_data.cur_view_labels = pred.views.labels;
            if gui_data.navigate_only_current_view_structure
                tmp = gui_data.cur_view_labels;
                gui_data.slices_to_navigate = find(tmp(gui_data.slice_id)==tmp);
            end
            tmp = gui_data.slices_to_navigate;
            gui_data.slices_to_navigate = circshift(tmp, -(find(tmp==gui_data.slice_id)-1));
        end

        while  ~gui_data.terminated && ~gui_data.recalculate && ~gui_data.force_plot_is_on
            gui_data.txt_handles(gui_data.status_txt_order_id).String =...
                "Waiting"+" - ("+round(fit_time,2)+" sec)";
            pause(0.1)
        end
        if gui_data.terminated
            return
        end
    end

    function graph_regen()
        gui_data.data_are_labeled = true;
        %======================== parse real datasets ===================
        %         [X,nodes_labels,views_labels] = parse_dataset(3);
        %         X = double(X);
        %         X = X(:,:,[3 6 2 8 9 1 4 5 7]);
        % %         [3 2 1 3 3 1 3 2 2]
        % %         views_labels = [1 1 2 2 2 3 3 3 3]
        %         for i=1:size(X,1)
        %             a(i)=sum(X(i,:,:),'all');
        %             i
        %         end
        %         X(a<=4,:,:)=[];
        %         X(:,a<=4,:)=[];
        %         size(X)
        % %         X = X(:,:,28);
        % %         views_labels = 0;
        % %         X=X(1:4:end,1:4:end,:);
        % % X= (X+permute(X,[2 1 3]))/2;
        %        cur_par = graph_tree_node;
        %        cur_par.loc = [1 1];
        %        cur_par.size = size(X,1);
        %===============================================================
        if ~exist('par_all','var');
            %%%%%%%%%%%%%%%%%%%%%%%%
%             par_all = graph_tree_root;
%             tmp = load('Datasets\openflights\tmp4.mat');
%             par_all.Data = tmp.X2;
%             tmp2 = string(tmp.airlines_final.region);
%             par_all.labels = tmp2; 
%             tmp2 = string(tmp.airports_final.region);
%             for i =1:numel(unique(tmp2))
%                 par_all.Children(i).labels = tmp2; 
%                 %
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%k
            %         
%             par_all = graph_tree_root; tmp = load('Datasets\reality_mining\extracted_data.mat');
%             par_all.Data = squeeze(sum(tmp.X,4));
%             par_all.labels = reshape(ones(size(par_all.Data,3),1).*[0:23],1,[]); 
%             par_all.Data = reshape(par_all.Data,size(par_all.Data,1),size(par_all.Data,2),[]);
%             for i =1:numel(unique(par_all.labels))
%                 if i == 1
%                     tmp2 = tmp.my_affil';
%                 else
%                     tmp2 = tmp.my_affil';
%                     %                             tmp2 = tmp.neighborhood';
%                 end
%                 tmp2 = strrep(string(tmp2),'_','\_');
%                 par_all.Children(i).labels = tmp2 
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%
            size_all= {...
                {30,20,10},
                {50,10},
                {10,50}
                };
            noise_level = 0.01;
%             sparsity_level_all = [0.94 0.93 0.92];
            sparsity_level_all = [1 1 1]*0.7;
            for k = 1:1%numel(sparsity_level_all)
                par_all = graph_tree_root;
                for i = 1:numel(size_all)
                    par_all.Children(i).is_symmetric = false;
                    par_all.Children(i).slices_num = 2;
                    par_all.Children(i).noise_level = noise_level;
                    par_all.Children(i).sparsity_level = sparsity_level_all(i);
                    for j = 1:numel(size_all{i})
                        par_all.Children(i).Children(j).type =...
                            'clique';
                        par_all.Children(i).Children(j).size = size_all{i}{j}*2;
                    end
                end
            end
        else
        end
        par_all.labels = reshape(par_all.labels,1,[]);
        for i=1:numel(par_all.Children)
            par_all.Children(i).labels = reshape(par_all.Children(i).labels,1,[]);
        end
        [X,par_all] = create_graph(par_all);
        tmp = unique(par_all.labels);
        par_all.labels = sum((tmp'==par_all.labels).*[1:numel(tmp)]',1);
        views_labels_names = "view label "+string(tmp);
        for i = 1:numel(par_all.Children)
            nodes_labels_names{i} = unique(par_all.Children(i).labels);
            par_all.Children(i).labels = sum((nodes_labels_names{i}'==par_all.Children(i).labels).*[1:numel(nodes_labels_names{i})]',1);
            nodes_labels_names{i} = "node label "+string(nodes_labels_names{i});
        end
        views_labels = par_all.labels;
        nodes_labels = {par_all.Children.labels};
        if ~gui_data.data_are_labeled
            views_labels = views_labels*0+1;
            for i = 1:numel(nodes_labels)
                nodes_labels{i}= nodes_labels{i}*0+1;
            end
        end

        L_rec = nan;

        if isempty(init_par)
            init_par = par_all.Children(1);
            cur_par = init_par;
        end

        % removes blank columns and rows
        %             tmp = 0;
        %             for i = 1:size(X,3)
        %                 XX = X(:,:,i);
        %                 tmp = tmp + sum(XX,1) + sum(XX,2)' - 2*diag(XX)';
        %             end
        %             tmp = find(tmp==0);
        %             X(tmp,:,:)=[];
        %             X(:,tmp,:)=[];
    end


    function  gui_plot()
        %turns off warning when plotting boxplots in tiledlayout
        warning('off','MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout')
        gui_data.community_types = {'clique','blank','star','triangular'};


        % Initializes gui. Runs only once
        if ~isfield(gui_data,'main_fig') ...
                || ~isprop(gui_data.main_fig,'Type') ...
                || ~strcmp(gui_data.main_fig.Type,'figure')
            gui_data.default_colors =     [
                0,       0.4470,  0.7410
                0.8500,  0.3250,  0.0980
                0.9290,  0.6940,  0.1250
                0.4940,  0.1840,  0.5560
                0.4660,  0.6740,  0.1880
                0.3010,  0.7450,  0.9330
                0.6350,  0.0780,  0.1840
                0,       0,       1
                0,       0.5,     0
                1,       0,       0
                0,       0.75,    0.75
                0.75,    0,       0.75
                0.75,    0.75,    0
                0.25,    0.25,    0.25];

            for i = 1:20
                gui_data.default_colors = [gui_data.default_colors;gui_data.default_colors.^5];
            end

            gui_data.init_X=X;

            gui_data.L_types_all=["Adjacency","Normalized Laplacian"];
            gui_data.clus_types_all = ["k-means","largest absolute value","large inner products"];

            gui_data.terminated = false;
            gui_data.recalculate = true;

            gui_data.plot_handles_is_on = true;
            gui_data.plot_orig_mat_is_on = true;
            gui_data.plot_orig_graph_is_on = true;
            gui_data.plot_clustered_mat_is_on = true;

            gui_data.init_X = X;
            gui_data.txt_handles = gobjects(0);
            gui_data.clustered_mat_inds = true(1,size(X,1));
            gui_data.cur_alg_opt_group_id = 7;
            gui_data.cur_alg_opt_id = 1;
            gui_data.slice_id = 1;
            gui_data.min_comp=1;
            gui_data.navigate_only_current_view_structure = true;
            gui_data.cur_view_labels = pred.views.labels;
            gui_data.slices_to_navigate = find(gui_data.cur_view_labels(gui_data.slice_id)==gui_data.cur_view_labels);
            gui_data.legend_visibility = "off";
            gui_data.navigate_based_on_orig_structure = false;

            info_bar_height = 0.1;


            gui_data.graph_modification_mode_is_on = false;
            gui_data.noise_modification_mode_is_on = false;
            gui_data.structure_modification_mode_is_on = false;
            gui_data.cluster_selection_mode_is_on = false;

            set(groot,  'defaultuicontrolbackgroundcolor',  [33,   39,   51]/255);
            set(groot,  'defaultuicontrolforegroundcolor',  [217,  215,  206]/255);
            set(groot,  'defaultuipanelbackgroundcolor',    [33,   39,   51]/255);
            set(groot,  'defaultfigurecolor',               [0     0     0])         
            set(groot,  'defaultaxescolor',                 [0.5   0.5   0.5]);
            set(groot,  'defaultaxesxcolor',                [1     1     1]);
            set(groot,  'defaultaxesycolor',                [1     1     1]);
            set(groot,  'defaultaxeszcolor',                [1     1     1]);
            set(groot,  'defaulttextcolor',                 [217,  215,  206]/255);
            set(groot,  'defaultaxesgridcolor',             [1     1     1]);
            set(groot,  'defaultUIControlFontSize',         10);
            gui_data.main_fig = figure('windowstate','maximized','MenuBar','None','keypressfcn',@key_press_callback,'units','normalized');
            gui_data.figures_panel = uipanel(gui_data.main_fig,'units','normalized','position',[0.0 0.00 1 1]);
            gui_data.axis_handles = gobjects(0);
            gui_data.settings_panel = uipanel(gui_data.main_fig,'units','normalized','position',[0.0 0.00 1 1]);
        end



        cur_X_slice = X(:,:,gui_data.slice_id);
        if gui_data.slice_id <= numel(pred.views.labels)
            k_calc = pred.views.labels(gui_data.slice_id);
            k_calc_ind = find(unique(pred.views.labels)==k_calc);
        else
            k_calc = nan;
        end


        k_orig_prev = k_orig;
        k_orig = views_labels(gui_data.slice_id);
        gui_data.k_orig_ind = find(unique(views_labels) == k_orig);
 


        if ~isnan(k_calc)
            tmp = perm_proj{k_calc_ind};
            cur_colors = tmp*gui_data.default_colors([1:size(tmp,2)],:);
        end


        if gui_data.plot_handles_is_on
            gui_data.alg_all={"","ComClus","Symmetric RichCom","CMNC"};
            gui_data.nodes_representation_types={'U','U*sqrt(B)','U*B'};

            GUI_txt = [];
            if  gui_data.graph_modification_mode_is_on
                if gui_data.noise_modification_mode_is_on == 1
                    GUI_txt{1} = "Graph Modification Mode: Noise";
                elseif gui_data.noise_modification_mode_is_on == 2
                    GUI_txt{1} = "Graph Modification Mode: Sparsity";
                elseif gui_data.structure_modification_mode_is_on
                    GUI_txt{1} = "Graph Modification Mode: Structure";
                else
                    GUI_txt{1} = "Graph Modification Mode: Type";
                end
            else
                GUI_txt{1} = "Graph Modification Mode: off";
            end

            i=1;
            while true
                if gui_data.cur_alg_opt_group_id ==i
                    tmp1 = ">"; tmp2 ="<";
                else
                    tmp1 =""; tmp2 ="";
                end
                switch i
                    case 1
                        GUI_txt{end+1} = "Normalization: " + tmp1+ gui_data.L_types_all(L_type_ind)+tmp2;
                    case 2
                        GUI_txt{end+1} = "Num. of Components: " + tmp1+R+tmp2;
                    case 3
                        GUI_txt{end+1} = "Clustering Method: " + tmp1+gui_data.clus_types_all(U_clus_type_ind)+tmp2;
                    case 4
                        if U_clus_type_ind==3
                            GUI_txt{end+1} = "Threshold: " + tmp1 + largest_inner_prod_thres + tmp2;
                        else    
                            if cl_num_auto
                                cl_num_U_aug = ["auto" string(U_cl_num)];
                            else
                                cl_num_U_aug = ["manual" string(U_cl_num)];
                            end
                            GUI_txt{end+1} = "Num. of Clusters: " + strjoin(cl_num_U_aug(1:min([gui_data.cur_alg_opt_id-1,end])),'  ') +"  "+tmp1+strjoin(cl_num_U_aug(gui_data.cur_alg_opt_id(gui_data.cur_alg_opt_id<=numel(cl_num_U_aug))),'  ')+tmp2+"  "+strjoin(cl_num_U_aug(gui_data.cur_alg_opt_id+1:end),'  ');
                        end
                    case 5
                        GUI_txt{end+1} = "Num. of Structures: " + tmp1+M+tmp2;
                    case 6
                        GUI_txt{end+1} = "Node Representation: "+ tmp1+gui_data.nodes_representation_types{U_nrm_type_ind(1)}+tmp2;
                    case 7
                        alg_opts_aug = [gui_data.alg_all{alg} string(alg_opts{alg})];
                        GUI_txt{end+1} = "Embedding Method: " + strjoin(alg_opts_aug(1:min([gui_data.cur_alg_opt_id-1,end])),'  ')+"  "+ tmp1+ ...
                            strjoin(alg_opts_aug(gui_data.cur_alg_opt_id(gui_data.cur_alg_opt_id<=numel(alg_opts_aug))),'  ') +tmp2+"  "+ ...
                            strjoin(alg_opts_aug(gui_data.cur_alg_opt_id+1:end),'  ');
                    case 8
                        GUI_txt{end+1} = "Thres: "+ tmp1+thres+tmp2;
                    case 9
                        GUI_txt{end+1} = "max_iters: "+ tmp1+max_iters+tmp2;
                    otherwise
                        gui_data.cur_alg_opt_group_id_max = i-1;
                        break;
                end
                i = i+1;
            end
            GUI_txt{end+1} = "Ploting";
            gui_data.status_txt_order_id = numel(GUI_txt);

            GUI_txt{end+1} = "Error:"+round(norm(L(:)-L_rec(:))/norm(L(:))*100,2)+" %" ;
            txt_autoresize(GUI_txt);
        end
        gui_data.txt_handles(gui_data.status_txt_order_id).String = "Ploting";

        ax_ind = 1;
        if gui_data.plot_orig_mat_is_on || gui_data.force_plot_is_on
            ax_ = get_axis(7);
            ax = get_axis(1);            
            if ~isfield(gui_data,'cur_com') || isempty(gui_data.cur_com)
                gui_data.cur_com = cur_par;
            end
            X_mask = ones(size(cur_X_slice));
            sz1 = gui_data.cur_com.size(1) + gui_data.cur_com.size(end)*gui_data.cur_com.is_bipartite;
            sz2 = gui_data.cur_com.size(end) + gui_data.cur_com.size(1)*gui_data.cur_com.is_bipartite;
            X_mask(gui_data.cur_com.loc(1)+gui_data.cur_com.offset+[1:sz1]-1,gui_data.cur_com.loc(2)+gui_data.cur_com.offset+[1:sz2]-1) = 0;
            XX = cur_X_slice- cur_X_slice.*X_mask/2;
            X_mask = ones(size(cur_X_slice,1),size(cur_X_slice,2),3);
            X_mask2 = zeros(size(cur_X_slice,1),size(cur_X_slice,2),3);
            inds_all =[];
            i = 1;
            for j = unique(nodes_labels{gui_data.k_orig_ind})
                inds_cur = sort(find(nodes_labels{gui_data.k_orig_ind}==j));
                X_mask(inds_cur,inds_cur,1) = gui_data.default_colors(i,1);
                X_mask(inds_cur,inds_cur,2) = gui_data.default_colors(i,2);
                X_mask(inds_cur,inds_cur,3) = gui_data.default_colors(i,3);
                X_mask2(inds_cur,inds_cur,:) = 1;
                inds_all = [inds_all inds_cur];
                i = i+1;
            end
            XX= XX.*X_mask;%+0.13*X_mask2.*(XX==0);
            tic
            if ~isfield(gui_data,'orig_mat_handle')
                gui_data.orig_mat_handle = imagesc(XX);
                ax.XTick = [];
                ax.YTick = [];
            else
                gui_data.orig_mat_handle.CData = XX;
            end
            if  gui_data.navigate_based_on_orig_structure && gui_data.navigate_only_current_view_structure
                tmp = ">";
            else
                tmp = "";
            end

            title(ax,"Original Adjacency Matrix - View "+gui_data.slice_id+"/"+size(X,3)+" -"+tmp+" Structure "+gui_data.k_orig_ind+"/"+numel(nodes_labels))
            axis(ax,'image')
            gui_data.plot_orig_mat_is_on = false;
            box(ax,'on')
            view_labels_image = ones(size(X,3),5,3);
            i = 1;
            for j = unique(views_labels)
                inds_cur = find(views_labels==j);
                view_labels_image(inds_cur,:,1) = gui_data.default_colors(i,1);
                view_labels_image(inds_cur,:,2) = gui_data.default_colors(i,2);
                view_labels_image(inds_cur,:,3) = gui_data.default_colors(i,3);
                i = i+1;
            end
            view_labels_image(gui_data.slice_id,3,:)= [0 0 0];
            imagesc(ax_,view_labels_image);
            ax_.XTick = [];
            ax_.YTick = [];
            gui_data.times_all(ax_ind)=toc;
%             disp("plot "+ax_ind+": "+gui_data.times_all(ax_ind))
            ax_ = get_axis(7);
            ax = get_axis(1);            
            ax.OuterPosition = [0.0 0.5 1/3-1/80 0.5];
            ax.Position = [0.0 0.5 1/3-1/80 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
            if exist("views_labels_names",'var') 
                add_legend(ax_,views_labels_names,0);
            end
            if exist("nodes_labels_names",'var') && iscell(nodes_labels_names) && gui_data.k_orig_ind<=numel(nodes_labels_names)
                tmp = gui_data.k_orig_ind;
                labels_names = nodes_labels_names{tmp}+"."+string(tmp);
                add_legend(ax,labels_names,0,"NorthEast");
            end
            ax_.OuterPosition = [1/3-1/80 0.5 1/80 0.5];
            ax_.Position = [1/3-1/80 0.5 1/80 0.5]-ax_.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
            tmp=0;
            if gui_data.legend_visibility=="on"
%                 ax_.Legend.Position(2)=1-ax_.Legend.Position(4);
%                 ax_.Legend.Position(1)=ax_.Position(1)-(ax_.Legend.Position(3)-2*ax_.Position(3))/2;
                tmp = ax_.Legend.Position(4);
            end
            ax_.OuterPosition = [1/3-1/80 0.5 1/80 max(0.5-tmp,0.5)];
            ax_.Position = [1/3-1/80 0.5 1/80 max(0.5-tmp,0.5)]-ax_.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
        end

        ax_ind = 4;
        % plot only nodes with non-zero degrees
        tmp = cur_X_slice~=0;
        tmp_inds = find((sum(tmp,1)>0)+(sum(tmp,2)>0)');
        %                     tmp_inds = 1:size(X,1);
        % normally should be put outside the "if recalulate", but is
        % put in here for efficiency. Plotting graphs is slow.
        if (gui_data.recalculate  && gui_data.plot_orig_graph_is_on) || gui_data.force_plot_is_on 
            tic

            ax = get_axis(4);
            %             ax = nexttile(gui_data.tl,2);
            gui_data.G_labeled = digraph(cur_X_slice(tmp_inds,tmp_inds),'omitselfloops');
            g_l = plot(ax,gui_data.G_labeled,'EdgeColor','k');
            title(ax,"Original Graph")
            i = 1;
            for j = unique(nodes_labels{gui_data.k_orig_ind})
                highlight(g_l,find(nodes_labels{gui_data.k_orig_ind}(tmp_inds)==j),'NodeColor',gui_data.default_colors(i,:))
                i = i+1;
            end
            gui_data.plot_orig_graph_is_on = false;
            gui_data.times_all(ax_ind)=toc;
%             disp("plot "+ax_ind+": "+gui_data.times_all(ax_ind))
            ax.OuterPosition = [0 0 1/3 0.5];
            ax.Position = [0 0 1/3 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
            if exist("nodes_labels_names",'var') && iscell(nodes_labels_names) && gui_data.k_orig_ind<=numel(nodes_labels_names)
                tmp = gui_data.k_orig_ind;
                labels_names = nodes_labels_names{tmp}+"."+string(tmp);
                add_legend(ax,labels_names,0);
            end
        end



        ax_ind = 3;
        if (gui_data.recalculate || gui_data.force_plot_is_on) && ~isnan(k_calc)
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ax = get_axis(3);
            %~~~~~~~~~~~~~~~~~~ plots first 3 dimensions of normalized latent representations ~~~~~~~~~~~~~~~~~
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %             for i=1:numel(U_nrm)
            %                 %         %Augment vectors with a third dimension if needed for plot3
            %                 %         P_aug = cat(2,P,zeros(size(U_nrm,1),3-size(U_nrm,2)));
            %                 %         pred_nodes_centroids_aug = cat(2,pred.nodes.centroids{i},zeros(size(pred.nodes.centroids{i},1),3-size(pred.nodes.centroids{i},2)));
            %                 U_nrm_aug{i} = cat(2,U_nrm{i},zeros(size(U_nrm{i},1),3-size(U_nrm{i},2)));
            %                 pred_nodes_centroids_nrm_aug{i} = cat(2,pred.nodes.centroids{i},zeros(size(pred.nodes.centroids{i},1),3-size(pred.nodes.centroids{i},2)));
            %             end
            %
            %             PP = U_nrm_aug{k_calc}(:,[gui_data.min_comp:gui_data.min_comp+2]);
            %             for j = unique(pred.nodes.labels{k_calc}(~isnan(pred.nodes.labels{k_calc})))'
            %                 zz = zeros(1,sum(pred.nodes.labels{k_calc}==j));
            %                 plot3([zz;PP(pred.nodes.labels{k_calc}==j,1)'],[zz;PP(pred.nodes.labels{k_calc}==j,2)'],[zz;PP(pred.nodes.labels{k_calc}==j,3)'],'Color',cur_colors(j,:));
            %                 hold on
            %             end
            %             grid on
            %             title("Latent representations - Comp.#: "+size(U_nrm{k_calc},2)+"/"+size(U,2))
            %             %     xlabel("Component " +  gui_data.min_comp)
            %             %     ylabel("Component " + (gui_data.min_comp+1))
            %             %     zlabel("Component " + (gui_data.min_comp+2))
            %             CC = pred_nodes_centroids_nrm_aug{k_calc}(:,[gui_data.min_comp:gui_data.min_comp+2])*0.5;
            %             zz = zeros(1,size(CC,1));
            %             plot3([zz;CC(:,1)'],[zz;CC(:,2)'],[zz;CC(:,3)'],'Color',[0 0 0],'LineWidth',3)
            %             axis equal
            %             hold off

            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %~~~~~~~~~~~~~~~ plots scatterplot of 2D MDS of normalized latent representations ~~~~~~~~~~~~~~~
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            U_nrm_tmp = U_nrm{k_calc};
            U_nrm_tmp(isnan(U_nrm_tmp))=0;
            P_MDS = cmdscale(squareform(pdist(U_nrm_tmp)),2);
            P_MDS_aug = cat(2,P_MDS,zeros(size(P_MDS,1),2-size(P_MDS,2)));
            node_IDs = 1:size(P_MDS,1);
            if k_calc <= numel(nodes_labels)
                tmp = nodes_labels{k_calc};
                tmp = (tmp'==unique(tmp));
                tmp = sum(tmp.*[1:size(tmp,2)],2)';
                scatter3(ax,P_MDS_aug(:,1),P_MDS_aug(:,2),node_IDs,[],gui_data.default_colors(tmp,:));
            else
                scatter3(ax,P_MDS_aug(:,1),P_MDS_aug(:,2),node_IDs,[]);
            end
            view(ax,2)
            axis(ax,'image')
            title(ax,"2D-MDS Latent representations - Comp.#: "+size(U_nrm{k_calc},2)+"/"+size(U,2))
            gui_data.times_all(ax_ind)=toc;
%             disp("plot "+ax_ind+": "+gui_data.times_all(ax_ind))
            ax.OuterPosition = [2/3 0.5 1/3 0.5];
            ax.Position = [2/3 0.5 1/3 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
            if exist("nodes_labels_names",'var') && iscell(nodes_labels_names) && gui_data.k_orig_ind<=numel(nodes_labels_names)
                tmp = k_calc_ind;
                if tmp<= numel(nodes_labels_names)
                    labels_names = nodes_labels_names{tmp}+"."+string(tmp);
                    add_legend(ax,labels_names,0,"NorthEast");
                end
            end
        end

        ax_ind = 2;
        ax_ = get_axis(8);
        if (gui_data.recalculate || gui_data.plot_clustered_mat_is_on || gui_data.force_plot_is_on) && ~isnan(k_calc)
            tic
            ax = get_axis(ax_ind);
            gui_data.inds_all=[];

            X_mask = ones(size(cur_X_slice));
            X_mask(gui_data.clustered_mat_inds,gui_data.clustered_mat_inds)=0;
            %             X_mask(gui_data.cur_com.loc(1)+gui_data.cur_com.offset+[1:gui_data.cur_com.size(1)]-1,gui_data.cur_com.loc(2)+gui_data.cur_com.offset+[1:gui_data.cur_com.size(end)]-1) = 0;
            XX = cur_X_slice-cur_X_slice.*X_mask/2;
            X_mask = ones(size(cur_X_slice,1),size(cur_X_slice,2),3);
            X_mask2 = zeros(size(cur_X_slice,1),size(cur_X_slice,2),3);


                X_mask3 = ones(max([1 round(size(cur_X_slice,2)/50)]),size(cur_X_slice,2),3);
            if k_calc<=numel(nodes_labels)
                i = 1;
                for j = unique(nodes_labels{k_calc})
                    inds_cur = find(nodes_labels{k_calc}==j);
                    X_mask3(:,inds_cur,1) = gui_data.default_colors(i,1);
                    X_mask3(:,inds_cur,2) = gui_data.default_colors(i,2);
                    X_mask3(:,inds_cur,3) = gui_data.default_colors(i,3);
                    i = i+1;
                end
            end

            gui_data.clustered_mat_cl_sizes = [];
            i = 1;
            for j = unique(pred.nodes.labels{k_calc_ind})'
                inds_cur = find(pred.nodes.labels{k_calc_ind}==j)';
                [~,tmp] = sort(nodes_labels{gui_data.k_orig_ind}(inds_cur));
                inds_cur = inds_cur(tmp);
                gui_data.clustered_mat_cl_sizes(end+1) = numel(inds_cur);
                X_mask(inds_cur,inds_cur,1) = cur_colors(i,1);
                X_mask(inds_cur,inds_cur,2) = cur_colors(i,2);
                X_mask(inds_cur,inds_cur,3) = cur_colors(i,3);
                X_mask2(inds_cur,inds_cur,:) = 1;
                gui_data.inds_all = [gui_data.inds_all inds_cur];
                i = i+1;
            end
            %             toc
            XX = X_mask.*XX+0.13*X_mask2.*(XX==0);
            XX2= XX(gui_data.inds_all,gui_data.inds_all,:);
            XX2 = cat(1,X_mask3(:,gui_data.inds_all,:),XX2);
            %             toc
            if isempty(clustered_matrix_handle)
                clustered_matrix_handle = imagesc(XX2);
                ax.XTick = [];
                ax.YTick = [];
                %             toc
                axis(ax,'image')
            else
                clustered_matrix_handle.CData = XX2;
            end

            %             toc
            if  ~gui_data.navigate_based_on_orig_structure && gui_data.navigate_only_current_view_structure
                tmp = ">";
            else
                tmp = "";
            end
            
            title(ax,"Clustered Adjacency Matrix -"+tmp+" Structure "+k_calc_ind+"/"+numel(unique(pred.views.labels)))
            
            if exist("nodes_labels_names",'var') && iscell(nodes_labels_names) && gui_data.k_orig_ind<=numel(nodes_labels_names)
                tmp = k_calc_ind;
                if tmp<= numel(nodes_labels_names)
                    labels_names = nodes_labels_names{tmp}+"."+string(tmp);
                    add_legend(ax,labels_names,0,"NorthEast");
                end
            end
            ax.OuterPosition = [1/3+1/80 0.5 1/3-1/80 0.5];
            ax.Position = [1/3+1/80 0.5 1/3-1/80 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
            gui_data.times_all(ax_ind)=toc;
            view_labels_image = ones(size(X,3),5,3);

            % Assumes the label ID of each predicted view cluster is the same
            % as its matching ground-truth view cluster
            for j = unique(pred.views.labels)
                inds_cur = find(pred.views.labels==j);
                view_labels_image(inds_cur,:,1) = gui_data.default_colors(j,1);
                view_labels_image(inds_cur,:,2) = gui_data.default_colors(j,2);
                view_labels_image(inds_cur,:,3) = gui_data.default_colors(j,3);
            end
            view_labels_image(gui_data.slice_id,3,:)= [0 0 0];
            imagesc(ax_,view_labels_image);
            ax_.XTick = [];
            ax_.YTick = [];
%             disp("plot "+ax_ind+": "+gui_data.times_all(ax_ind))
            gui_data.plot_clustered_mat_is_on = false;
        end
            ax_7 = get_axis(7);
        ax_.OuterPosition = [1/3 0.5 1/80 ax_7.Position(4)];
        ax_.Position = [1/3 0.5 1/80  ax_7.Position(4)]-ax_.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];

        ax_ind = 5;
        ax = get_axis(5);
        if (gui_data.recalculate || gui_data.force_plot_is_on) && ~isnan(k_calc)
            tic
            gui_data.G = digraph(cur_X_slice(tmp_inds,tmp_inds),'omitselfloops');
            g1 = plot(ax,gui_data.G,'EdgeColor','k','NodeColor','k');
            title(ax,"Clustered Graph")
            i = 1;
            for j = unique(pred.nodes.labels{k_calc_ind})'
                highlight(g1,find(pred.nodes.labels{k_calc_ind}(tmp_inds)==j),'NodeColor',cur_colors(i,:));
                i = i+1;
            end
            gui_data.times_all(ax_ind)=toc;
%             disp("plot "+ax_ind+": "+gui_data.times_all(ax_ind))
            ax.OuterPosition = [1/3 0 1/3 0.5];
            ax.Position = [1/3 0 1/3 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1]; 
            if exist("nodes_labels_names",'var') && iscell(nodes_labels_names) && gui_data.k_orig_ind<=numel(nodes_labels_names) 
                tmp = k_calc_ind;
                if tmp<= numel(nodes_labels_names)
                    labels_names = nodes_labels_names{tmp}+"."+string(tmp);
                    add_legend(ax,labels_names,0);
                end
            end
        end
        ax_ind = 6;
        if ( gui_data.recalculate || gui_data.force_plot_is_on ) && false
            ax = get_axis(6);
            tic
            %             plot(thres_all,nodes_scan_qual,"o-")
            for i = 1:size(nodes_scan_qual,3)
                prctile_plot(ax,thres_all,nodes_scan_qual(:,:,i));
                hold(ax,'on')
            end
            prctile_plot(ax,thres_all,mean(nodes_scan_qual,3),'Color','black');
            prctile_plot(ax,thres_all,views_scan_qual,'Color','black','LineStyle','--');
            hold(ax,'off')
            grid(ax,'on')
            xlabel(ax,"Inner product threshold")
            title(ax,"Clustering quality")
            ylabel(ax,selected_param_scan_name)
            gui_data.times_all(ax_ind)=toc;
            ylim(ax,[-1.1 1.1]);
            %             legend(ax,"View cluster "+[1:size(nodes_scan_qual,3)],'Location','BestOutside')

            labels_names = ["views","nodes - average", "nodes - "+views_labels_names]; 
            add_legend(ax,labels_names,1);
%             disp("plot "+ax_ind+": "+gui_data.times_all(ax_ind))
            %             if ~isempty(eigvals)
            %                 plot(eigvls,"o-")
            %                 title("Eigenalues of all estimated structures")
            %                 ylim([0 1])
            %             else
            %                 ax.delete;
            %             end

            ax.OuterPosition = [2/3 0 1/3 0.5];
            ax.Position = [2/3 0 1/3 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
        end
%         disp("All plots: "+sum(gui_data.times_all))
        gui_data.main_fig.SizeChangedFcn=@GUI_auto_resize; 

        if gui_data.recalculate || gui_data.force_plot_is_on
            % plots U and V separately
            %     ax_ind = ax_ind+1;
            %     nexttile(ax_ind)
            %     zz = zeros(1,size(V_t,1));
            %     U_t_aug = cat(2,U_t,zeros(size(U_t,1),3-size(U_t,2)));
            %     V_t_aug = cat(2,V_t,zeros(size(V_t,1),3-size(V_t,2)));
            %     plot3([zz;U_t_aug(:,1)'],[zz;U_t_aug(:,2)'],[zz;U_t_aug(:,3)'],'Color',cur_colors(1,:));
            %     hold on
            %     plot3([zz;V_t_aug(:,1)'],[zz;V_t_aug(:,2)'],[zz;V_t_aug(:,3)'],'Color',cur_colors(2,:));
            %     hold off
            %     axis equal

            %
            %
            %             if ~isempty(eigvals)
            %                 ax_ind = ax_ind+1;
            %                 nexttile(ax_ind)
            %                 %using plots - FAST
            %                 data = eigvals;
            %                 %             if ~isfield(gui_data,'st1_handle') || size(data,2)~=numel(gui_data.st1_handle)
            %
            %                 gui_data.st1_handle = plot(sort(eig(L),'descend'),'-xk');
            %                 hold on
            %                 for i = 1:size(data,2)
            %                     gui_data.st1_handle = plot([1:size(data,1)]+(0.2)*i,data(:,i),'-o','Color',cur_colors(i,:));
            %                 end
            %                 xlim([1 30])
            %                 hold off;
            %                 grid on
            %                 title("Eigenvalues of clusters")
            %                 %             else
            %                 %
            %                 %                         end
            %
            %                 %using stackedplots - SLOW
            %                 %             data = eigvals(1:min([size(eigvals,1) 5]),:).^2;
            %                 %             if ~isfield(gui_data,'st1_handle') || size(data,2)~=size(gui_data.st1_handle.YData,2)
            %                 %                 gui_data.st1_handle = stackedplot(data);
            %                 %             else
            %                 %                 gui_data.st1_handle.YData(1:size(data,1),:) = data;
            %                 %             end
            %                 %
            %                 %
            %                 %             tic
            %                 %             tmp_lims = cell2mat({gui_data.st1_handle.AxesProperties.YLimits}');
            %                 %
            %                 %             for i = 1:numel(gui_data.st1_handle.AxesProperties)
            %                 %                 gui_data.st1_handle.AxesProperties(i).YLimits(1)= min(tmp_lims(:,1));
            %                 %                 gui_data.st1_handle.AxesProperties(i).YLimits(2)= max(tmp_lims(:,2));
            %                 %             end
            %                 %             toc
            %
            %
            %                 gui_data.times_all(ax_ind)=toc;
            %                 disp("plot "+ax_ind+": "+gui_data.times_all(ax_ind))
            %             end

            %             tic
            %             ax_ind = ax_ind+1;
            %             nexttile(ax_ind)
            %
            %             %commented out because it is buggy
            %             %             if ~isfield(gui_data,'st2_handle') || size(U_nrm,2)~=size(gui_data.st2_handle.YData,2)
            %             gui_data.st2_handle = stackedplot(U_nrm);
            %             %             else
            %             %                 gui_data.st2_handle.YData = U_nrm;
            %             %             end
            %
            %             gui_data.times_all(ax_ind)=toc;
            %             disp("plot "+ax_ind+": "+gui_data.times_all(ax_ind))
            %

            %
            %
            %         ax_ind = ax_ind+1;
            %         ax(ax_ind)= nexttile(gui_data.tl,ax_ind);
            %         boxplot(cl_purity_1)
            %         h = flip(findobj(ax(ax_ind).Children(1).Children,'Tag','Box'));
            %         for j=1:length(h)
            %             patch(get(h(j),'XData'),get(h(j),'YData'),cur_colors(j,:));
            %         end
            %         ylim([0 1])
            %         set(ax(ax_ind),'children',flipud(get(ax(ax_ind),'children')))
            %
            %
            %         ax_ind = ax_ind+1;
            %         ax(ax_ind)= nexttile(ax_ind);
            %         boxplot(cl_purity_2)
            %         h = flip(findobj(ax(ax_ind).Children(1).Children,'Tag','Box'));
            %         for j=1:length(h)
            %             patch(get(h(j),'XData'),get(h(j),'YData'),cur_colors(j,:));
            %         end
            %         ylim([0 1])
            %         set(ax(ax_ind),'children',flipud(get(ax(ax_ind),'children')))
            %
        end
        GUI_auto_resize;
    end

    function leg = add_legend(ax,labels_names,type,loc)
        if gui_data.legend_visibility=="on"
            if ~exist('loc','var')
                loc = 'Best';
            end
            hold(ax,'on')
            legend_dummy_plots = [];
            if type==1
                legend_dummy_plots = num2cell(ax.Children);
            else
                for i = 1:numel(labels_names)
                    legend_dummy_plots{i} = plot(ax,nan,'.','MarkerSize',70,'Color',gui_data.default_colors(i,:));
                end
            end
            hold(ax,'off')
            leg = legend(ax,[legend_dummy_plots{:}],cellstr(labels_names)','FontSize',10,'Location',loc,'AutoUpdate','off');
            setappdata(ax,'LegendColorbarManualSpace',1);
            setappdata(ax,'LegendColorbarReclaimSpace',1);
            set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.2;.2;.2;.8]))
        end
    end

    function txt_autoresize(GUI_txt)
        tmp =0 ;
        tmp3 = 0;
        for i = 1 : numel(GUI_txt)
            tmp3 = tmp3+strlength(GUI_txt{i});                % fix condition
        end
        tmp0 =  numel(gui_data.txt_handles)<numel(GUI_txt);
        max_string_height = 0;
        for i = 1 : numel(GUI_txt)
            %                     gui_data.txt_handles(i) = uicontrol('Parent',gui_data.settings_panel,'String',GUI_txt{i},'Style', 'text','units','normalized','outerposition',[1/numel(GUI_txt)*(i-1) 0  1/numel(GUI_txt)*0.9 0.8 ]);
            if tmp0
                gui_data.txt_handles(i) = uicontrol('Parent',gui_data.settings_panel,'String',GUI_txt{i},'Style', 'text','Units','normalized','Position',[0 0 1 1]);
            end
            gui_data.txt_handles(i).String = GUI_txt{i};
            gui_data.txt_handles(i).Units = 'normalized';
            if tmp0
                gui_data.txt_handles(i).Position(1) = tmp;
                gui_data.txt_handles(i).Position(3) = max([20,strlength(GUI_txt{i})])/tmp3;
                tmp = tmp + gui_data.txt_handles(i).Position(3); 
                gui_data.txt_handles(i).Position(2) = 0;
            end
        end
        tmp4 = sum(gui_data.txt_handles(end).Position([1,3]));
        for i = 1 : numel(gui_data.txt_handles)
            cur_txt_handle = gui_data.txt_handles(i);
            cur_txt_handle.Units = 'characters';
            tmp2 = cur_txt_handle.Extent(4);
            cur_txt_handle.Units = 'pixels';
            char_height = cur_txt_handle.Extent(4)/tmp2;
            max_string_height = max([max_string_height, (ceil(cur_txt_handle.Extent(3)/cur_txt_handle.Position(3))+1)*char_height]);
            cur_txt_handle.Units = 'normalized';
            cur_txt_handle.Position([1,3]) = cur_txt_handle.Position([1,3])/tmp4;
        end
        gui_data.settings_panel.Units = 'pixels';
        gui_data.settings_panel.Position(4) = max_string_height;
        gui_data.settings_panel.Units = 'normalized';
        tmp5 = gui_data.settings_panel.Position(4);
        gui_data.figures_panel.Position([2,4]) = [tmp5,1-tmp5];
        gui_data.plot_handles_is_on = false;
    end

    function GUI_auto_resize(src,evnt)
        drawnow
        ax_ = get_axis(7);
        ax = get_axis(1);            
        ax.OuterPosition = [0.0 0.5 1/3-1/80 0.5];
        ax.Position = [0.0 0.5 1/3-1/80 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];

        tmp=0;
        if gui_data.legend_visibility=="on"
            ax_.Legend.Position(2)=1-ax_.Legend.Position(4);
            ax_.Legend.Position(1)=ax_.Position(1)-(ax_.Legend.Position(3)-2*ax_.Position(3))/2;
            tmp = ax_.Legend.Position(4);
        end
        ax_.OuterPosition = [1/3-1/80 0.5 1/80 0.5-tmp*0];
        ax_.Position = [1/3-1/80 0.5 1/80 0.5-tmp*0]-ax_.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];

        ax = get_axis(4);
        ax.OuterPosition = [0 0 1/3 0.5];
        ax.Position = [0 0 1/3 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
        ax = get_axis(3);
        ax.OuterPosition = [2/3 0.5 1/3 0.5];
        ax.Position = [2/3 0.5 1/3 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
        ax_ = get_axis(8);
        ax = get_axis(2);
        ax.OuterPosition = [1/3+1/80 0.5 1/3-1/80 0.5];
        ax.Position = [1/3+1/80 0.5 1/3-1/80 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
        ax_7 = get_axis(7);
        ax_.OuterPosition = [1/3 0.5 1/80 ax_7.Position(4)];
        ax_.Position = [1/3 0.5 1/80  ax_7.Position(4)]-ax_.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
        ax = get_axis(5);
        ax.OuterPosition = [1/3 0 1/3 0.5];
        ax.Position = [1/3 0 1/3 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
        ax = get_axis(6);
        ax.OuterPosition = [2/3 0 1/3 0.5];
        ax.Position = [2/3 0 1/3 0.5]-ax.TightInset*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1];
        txt_autoresize([]);
    end

    function update_slices_to_navigate()
        if ~gui_data.navigate_only_current_view_structure
            gui_data.slices_to_navigate = 1:size(X,3);
        else
            gui_data.slices_to_navigate = find(gui_data.cur_view_labels(gui_data.slice_id)==gui_data.cur_view_labels)
        end
        gui_data.slices_to_navigate = circshift(gui_data.slices_to_navigate, -(find(gui_data.slices_to_navigate==gui_data.slice_id)-1));
    end

    function  next_item = get_next_item(cur_item,all_items)
        next_item =  all_items(circshift(cur_item == all_items,1));
    end

    function  prev_item = get_prev_item(cur_item,all_items)
        prev_item =  all_items(circshift(cur_item == all_items,-1));
    end

    function key_press_callback(h,evnt)
        graph_regen = @graph_regen;
        update_slices_to_navigate = @update_slices_to_navigate;
        gui_plot = @gui_plot;
        get_next_item = @get_next_item;
        get_prev_item = @get_prev_item;
        [tmp,tmp1,tmp2,tmp3,tmp4,tmp5,cur_com_type_ind,cur_ind] = deal([]);
        run key_press_callback.m
    end
    
    function ax = get_axis(ax_ind)
        if numel(gui_data.axis_handles)>=ax_ind && isgraphics(gui_data.axis_handles(ax_ind),'axes')
            ax = gui_data.axis_handles(ax_ind);
        else
            ax = axes(gui_data.figures_panel);
            gui_data.axis_handles(ax_ind) = ax;
        end
    end
end
