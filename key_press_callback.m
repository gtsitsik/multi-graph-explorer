% Definition of function single_run_GUI/key_press_callback
switch evnt.Key
    case 'alt'
        keyboard;
    case 'f1'
        gui_data.main_fig.MenuBar = get_next_item(gui_data.main_fig.MenuBar,["none","figure"]);
    case 'space'
        gui_data.recalculate = true;
        gui_data.cluster_selection_mode_is_on = false;
        gui_data.plot_handles_is_on = true;
        return
    case 'g'
        if ~gui_data.graph_modification_mode_is_on
            gui_data.graph_modification_mode_is_on = true;
        else
            gui_data.graph_modification_mode_is_on = false;
            gui_data.noise_modification_mode_is_on = false;
            gui_data.structure_modification_mode_is_on = false;
        end
        gui_data.cur_com = cur_par;
        gui_data.plot_orig_mat_is_on = true;
    case 'h'
        gui_data.legend_visibility = get_next_item(gui_data.legend_visibility,["on","off"]);
        if gui_data.legend_visibility == "off"
            delete(findobj('Type','legend'))
        end
        gui_data.force_plot_is_on = true;
    case {'quote','tab'}
        % TODO: solve conflict with other use of tab
        if evnt.Key == "tab"
            gui_data.navigate_based_on_orig_structure =~gui_data.navigate_based_on_orig_structure;
            if gui_data.navigate_based_on_orig_structure
                gui_data.cur_view_labels = views_labels;
            else
                gui_data.cur_view_labels = pred.views.labels;
            end
        end
        if evnt.Key == "quote"
            gui_data.navigate_only_current_view_structure = ~gui_data.navigate_only_current_view_structure;
        end
        update_slices_to_navigate();
        gui_data.plot_orig_mat_is_on = true;
        gui_data.plot_clustered_mat_is_on = true;
    case {'comma','period','slash'}
        switch evnt.Key
            case {'comma','period'}

                switch evnt.Key
                    case 'comma'
                        gui_data.slices_to_navigate = circshift(gui_data.slices_to_navigate,1);
                    case 'period'
                        gui_data.slices_to_navigate = circshift(gui_data.slices_to_navigate,-1);
                end
                gui_data.slice_id = gui_data.slices_to_navigate(1);
            case 'slash'
                tmp = gui_data.cur_view_labels;
                tmp2 = unique(tmp);
                gui_data.slice_id = find(tmp2(circshift(tmp(gui_data.slice_id)==tmp2,1))==tmp,1);

                if gui_data.navigate_only_current_view_structure
                    gui_data.slices_to_navigate = find(tmp(gui_data.slice_id)==tmp);
                end
        end
        % TODO: assumes that the view label of view structure i is a larger number than that of view structure j when i>j.
        if gui_data.data_are_labeled
            cur_par = par_all.Children(find(views_labels(gui_data.slice_id)==unique(views_labels)));
            gui_data.cur_com = cur_par;
        end
        gui_data.plot_orig_mat_is_on = true;
        gui_data.plot_orig_graph_is_on = true;
        gui_data.plot_clustered_mat_is_on = true;
        gui_data.plot_clustered_graph_is_on = true;

    case 'escape'
        gui_data.terminated = true;
        close(gui_data.main_fig)
        reset(groot)
        return
end
if ~gui_data.graph_modification_mode_is_on
    switch evnt.Key
        case 'return'
            gui_data.force_plot_is_on = true;
        case 'rightarrow'
            R = R+1;
        case 'leftarrow'
            R = max([R-1 0]);
        case {'add','subtract'}
            switch gui_data.cur_alg_opt_group_id
                case 4
                    switch U_clus_type_ind
                        case 3
                            switch evnt.Key
                                case 'add'
                                    largest_inner_prod_step = min([largest_inner_prod_step*10 0.1]);
                                case 'subtract'
                                    largest_inner_prod_step = max([largest_inner_prod_step/10 1e-5]);
                            end
                    end
                case 5
                    switch evnt.Key
                        case 'add'  
                            M = M+1;
                            U_cl_num(end+1) =1;
                            %                     R = sum(U_cl_num);
                        case 'subtract'
                            if M-1>0
                                M = M-1;
                                U_cl_num(end) =[];
                                %     R = sum(U_cl_num);
                            end
                    end
            end
        case 'l'
            L_type_ind = mod(L_type_ind,numel(gui_data.L_types_all))+1;
            %                     gui_data.recalculate = true;
        case '7'
            gui_data.min_comp = mod(gui_data.min_comp,(R-2)*(R>2)+(R<=2))+1;
        case '9'
            gui_data.min_comp = max([gui_data.min_comp-1 1]);
        case {'q','e','d','a'}

            switch gui_data.cur_alg_opt_group_id
                case {1,2,3,5,6,8,9}
                    gui_data.cur_alg_opt_id_max = 1;
                case 4
                    gui_data.cur_alg_opt_id_max = numel(U_cl_num)+1;
                case 7
                    gui_data.cur_alg_opt_id_max = numel(alg_opts{alg})+1;
            end
            switch evnt.Key
                case 'q'
                    gui_data.cur_alg_opt_id = 1;
                    gui_data.cur_alg_opt_group_id = gui_data.cur_alg_opt_group_id-1;
                    if gui_data.cur_alg_opt_group_id == 0
                        gui_data.cur_alg_opt_group_id = gui_data.cur_alg_opt_group_id_max;
                    end
                case 'e'
                    gui_data.cur_alg_opt_id = 1;
                    gui_data.cur_alg_opt_group_id = mod(gui_data.cur_alg_opt_group_id,gui_data.cur_alg_opt_group_id_max)+1;
                case 'd'
                    gui_data.cur_alg_opt_id = mod(gui_data.cur_alg_opt_id,gui_data.cur_alg_opt_id_max)+1;
                case 'a'
                    gui_data.cur_alg_opt_id = gui_data.cur_alg_opt_id-1;
                    if gui_data.cur_alg_opt_id  == 0
                        gui_data.cur_alg_opt_id = gui_data.cur_alg_opt_id_max;
                    end
            end


        case 'w'
            switch gui_data.cur_alg_opt_group_id
                case 1
                    L_type_ind = mod(L_type_ind,numel(gui_data.L_types_all))+1;
                case 2
                    R = R+1;
                case 3
                    U_clus_type_ind = mod(U_clus_type_ind,numel(gui_data.clus_types_all))+1;
                case 4
                    switch U_clus_type_ind 
                        case 3
                            largest_inner_prod_thres = min([largest_inner_prod_thres+largest_inner_prod_step,1]);
                        otherwise
                            switch gui_data.cur_alg_opt_id
                                case 1
                                    cl_num_auto=~cl_num_auto;
                                otherwise
                                    % use to change U_cl_num corresponding to the structure
                                    % of the current slice, given labeled data
                                    U_cl_num(gui_data.cur_alg_opt_id-1) = min([ (U_cl_num(gui_data.cur_alg_opt_id-1)+1), size(gui_data.default_colors,1)]);
                            end
                    end
                case 5
                    M = M+1;
                    U_cl_num(end+1) =1;
                    %                             R = sum(U_cl_num);
                case 6
                    U_nrm_type_ind(1)=mod(U_nrm_type_ind(1),numel(gui_data.nodes_representation_types))+1;
                case 7
                    switch gui_data.cur_alg_opt_id
                        case 1
                            alg = mod(alg,numel(gui_data.alg_all))+1;
                            if alg==1
                                alg = mod(alg,numel(gui_data.alg_all))+1;
                            end
                            gui_data.cur_alg_opt_id_max = numel(alg_opts{alg})+1;
                        otherwise
                            switch alg
                                case 1
                                    switch gui_data.cur_alg_opt_id
                                        case 4
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1} = alg_opts{alg}{gui_data.cur_alg_opt_id-1}+0.1;
                                        case {2,3}
                                            switch gui_data.cur_alg_opt_id
                                                case 2
                                                    tmp2 = 'U+1R';
                                                case 3
                                                    tmp2 = 'U+1';
                                            end
                                            tmp = alg_opts{alg}{gui_data.cur_alg_opt_id-1};
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1} = tmp2(circshift(tmp==tmp2,-1));
                                    end
                                case 2
                                    switch gui_data.cur_alg_opt_id
                                        case 4
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1} = alg_opts{alg}{gui_data.cur_alg_opt_id-1}*10;
                                        otherwise %ComClus rho
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1} = alg_opts{alg}{gui_data.cur_alg_opt_id-1} + 0.01;
                                    end
                                case {3,4} %Richcom rho | CMNC delta
                                    switch gui_data.cur_alg_opt_id
                                        case 2
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1} = alg_opts{alg}{gui_data.cur_alg_opt_id-1} + 0.01;
                                        case 3

                                            tmp2 = ["random","true"];

                                            tmp = alg_opts{alg}{gui_data.cur_alg_opt_id-1};
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1} = tmp2(circshift(tmp==tmp2,-1));
                                    end
                            end
                    end
                case 8
                    thres = thres*10;
                case 9
                    max_iters = max_iters+100;
            end
            %                     gui_data.recalculate = true;
        case 's'
            switch gui_data.cur_alg_opt_group_id
                case 1
                    L_type_ind = L_type_ind-1;
                    if L_type_ind ==0
                        L_type_ind = numel(gui_data.L_types_all);
                    end
                case 2
                    R = max([R-1 1]);
                case 3
                    U_clus_type_ind = U_clus_type_ind-1;
                    if U_clus_type_ind ==0
                        U_clus_type_ind = numel(gui_data.clus_types_all);
                    end
                case 4
                    switch U_clus_type_ind 
                        case 3
                            largest_inner_prod_thres = max([largest_inner_prod_thres-largest_inner_prod_step,0]);
                        otherwise
                            switch gui_data.cur_alg_opt_id
                                case 1
                                    cl_num_auto=~cl_num_auto;
                                otherwise
                                    % use to change U_cl_num corresponding to the structure
                                    % of the current slice, given labeled data
                                    U_cl_num(gui_data.cur_alg_opt_id-1) = max([U_cl_num(gui_data.cur_alg_opt_id-1)-1 1]);
                            end
                    end
                case 5
                    if M-1>0
                        M = M-1;
                        U_cl_num(end) =[];
                        %                                 R = sum(U_cl_num);
                    end
                case 6
                    U_nrm_type_ind(1)= U_nrm_type_ind(1)-1;
                    if U_nrm_type_ind(1) ==0
                        U_nrm_type_ind(1) = numel(gui_data.nodes_representation_types);
                    end
                case 7
                    switch gui_data.cur_alg_opt_id
                        case 1
                            alg = alg-1;
                            if alg  == 1
                                alg = numel(gui_data.alg_all);
                            end
                            gui_data.cur_alg_opt_id_max = numel(alg_opts{alg})+1;
                        otherwise
                            switch alg
                                case 1
                                    switch gui_data.cur_alg_opt_id
                                        case 4
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1}= max([alg_opts{alg}{gui_data.cur_alg_opt_id-1}-0.1,0]);
                                        case {2,3}
                                            switch gui_data.cur_alg_opt_id
                                                case 2
                                                    tmp2 = 'U+1R';
                                                case 3
                                                    tmp2 = 'U+1';
                                            end
                                            tmp = alg_opts{alg}{gui_data.cur_alg_opt_id-1};
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1} = tmp2(circshift(tmp==tmp2,1));
                                    end
                                case 2
                                    switch gui_data.cur_alg_opt_id
                                        case 4
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1} = alg_opts{alg}{gui_data.cur_alg_opt_id-1}/10;
                                        otherwise %ComClus rho
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1} = max([alg_opts{alg}{gui_data.cur_alg_opt_id-1}-0.01 , 0]);
                                    end
                                case {3,4} %Richcom rho | CMNC delta
                                    switch gui_data.cur_alg_opt_id
                                        case 2
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1} = max([alg_opts{alg}{gui_data.cur_alg_opt_id-1}-0.01 , 0]);
                                        case 3
                                            tmp2 = ["random","true"];
                                            tmp = alg_opts{alg}{gui_data.cur_alg_opt_id-1};
                                            alg_opts{alg}{gui_data.cur_alg_opt_id-1} = tmp2(circshift(tmp==tmp2,1));
                                    end
                            end
                    end
                case 8
                    thres = thres/10;
                case 9
                    max_iters = max_iters-100;
            end
        case 'u'
            U_nrm_type_ind(1)=mod(U_nrm_type_ind(1),3)+1;
            %                     gui_data.recalculate = true;
    end
else
    if gui_data.cluster_selection_mode_is_on
        switch evnt.Key
            case {'tab','rightarrow','leftarrow'}
                switch evnt.Key
                    case 'tab'
                        gui_data.cluster_selection_mode_is_on = false;
                        gui_data.clustered_mat_inds = gui_data.inds_all;
                    otherwise
                        switch evnt.Key
                            case 'rightarrow'
                                gui_data.clustered_mat_cur_cluster = mod(gui_data.clustered_mat_cur_cluster,numel(gui_data.clustered_mat_cl_sizes))+1;
                            case 'leftarrow'
                                gui_data.clustered_mat_cur_cluster = gui_data.clustered_mat_cur_cluster-1;
                                if  gui_data.clustered_mat_cur_cluster==0
                                    gui_data.clustered_mat_cur_cluster = numel(gui_data.clustered_mat_cl_sizes);
                                end
                        end
                        tmp1 = gui_data.clustered_mat_cl_sizes;
                        tmp2 = gui_data.clustered_mat_cur_cluster;
                        gui_data.clustered_mat_inds = gui_data.inds_all([1:tmp1(tmp2)]+sum(tmp1(1:tmp2-1)));
                end
                gui_data.plot_clustered_mat_is_on = true;
            case 'c'
                %                         cur_par = gui_data.cur_com;
                %                         graph_regen();
                % X= (gui_data.init_X(1:100,1:100));
                X_mask = zeros(size(gui_data.init_X));
                X_mask(gui_data.clustered_mat_inds,gui_data.clustered_mat_inds) = 1;
                X=gui_data.init_X.*X_mask;
                gui_data.plot_orig_mat_is_on = true;
                gui_data.plot_orig_graph_is_on = true;
        end
    else
        gui_data.plot_orig_mat_is_on = true;
        if gui_data.noise_modification_mode_is_on
            switch gui_data.noise_modification_mode_is_on
                case 1
                    tmp = gui_data.cur_com.noise_level;
                case 2
                    tmp = gui_data.cur_com.sparsity_level;
            end
            switch evnt.Key
                case {'add','subtract','d','0','1','2','3','4','5','6','7','8','9','numpad0','numpad1','numpad2','numpad3','numpad4','numpad5','numpad6','numpad7','numpad8','numpad9'}
                    switch evnt.Key
                        case 'add'
                            if tmp>0
                                if tmp+0.01<=1
                                    tmp = tmp+0.01;
                                end
                            else
                                tmp=0.01;
                            end
                        case 'subtract'
                            tmp = max([tmp-0.01 0]);
                        case 'd'
                            if gui_data.noise_modification_mode_is_on == 1
                                gui_data.cur_com.has_global_noise = ~gui_data.cur_com.has_global_noise;
                            elseif gui_data.noise_modification_mode_is_on == 2
                                gui_data.cur_com.has_global_sparsity = ~gui_data.cur_com.has_global_sparsity;
                            end

                        case {'0','1','2','3','4','5','6','7','8','9','numpad0','numpad1','numpad2','numpad3','numpad4','numpad5','numpad6','numpad7','numpad8','numpad9'}
                            evnt.Key(end)
                            tmp = str2double(evnt.Key(end))/10;
                    end
                    switch gui_data.noise_modification_mode_is_on
                        case 1
                            gui_data.cur_com.noise_level = tmp;
                        case 2
                            gui_data.cur_com.sparsity_level = tmp;
                    end
                    graph_regen();
                    gui_data.init_X=X;
                    gui_data.plot_orig_graph_is_on = true;
            end
        elseif gui_data.structure_modification_mode_is_on
            switch evnt.Key
                case 'i'
                    cur_ind = gui_data.k_orig_ind;
                    par_all.Children(cur_ind+2:end+1) = ...
                        par_all.Children(cur_ind+1:end)
                    par_all.Children(cur_ind+1) = graph_tree_node;
                    par_all.Children(cur_ind+1).noise_level = 0.001;
                    par_all.Children(cur_ind+1).is_symmetric = false;
                    par_all.Children(cur_ind+1).Children(1).size = [100];
                    par_all.Children(cur_ind+1).Children(1).type = 'clique';
                    par_all.Children(cur_ind+1).Children(1).size = [10];
                    par_all.Children(cur_ind+1).slices_num = 2;
                    graph_regen();
                    update_slices_to_navigate();
                    gui_data.init_X=X;
                    gui_data.plot_orig_graph_is_on = true;
                case 'd'
                    % TODO: Verify correctness
                    tmp = views_labels(gui_data.slice_id);
                    tmp2 = unique(views_labels);
                    tmp3 = tmp==tmp2;
                    tmp4 = circshift(tmp3,-1);
                    tmp5 = tmp2(tmp4);
                    gui_data.slice_id = 1; 
                    if gui_data.data_are_labeled
                        cur_par = par_all.Children(tmp4);
                        gui_data.cur_com = cur_par;
                    end
                    par_all.Children(tmp3)=[];
                    views_labels(views_labels==tmp)=[];
                    graph_regen();
                    update_slices_to_navigate();
                case 'add'
                    cur_par.slices_num = cur_par.slices_num+1;
                    graph_regen();
                    update_slices_to_navigate();
                case 'subtract'
                    if cur_par.slices_num-1>0
                        cur_par.slices_num = cur_par.slices_num-1;
                        if gui_data.slice_id == size(X,3)
                            gui_data.slice_id = gui_data.slice_id-1;
                        end
                        graph_regen();
                        update_slices_to_navigate();
                    end
            end
        else
            switch evnt.Key
                case 'uparrow'
                    cur_com_type_ind = find(strcmp(gui_data.community_types,gui_data.cur_com.type));
                    if ~isempty(cur_com_type_ind)
                        gui_data.cur_com.type = gui_data.community_types{mod(cur_com_type_ind,numel(gui_data.community_types))+1};
                        graph_regen();
                        gui_data.init_X=X;
                        gui_data.plot_orig_graph_is_on = true;
                    end
                case 'downarrow'
                    cur_com_type_ind = find(strcmp(gui_data.community_types,gui_data.cur_com.type));
                    if ~isempty(cur_com_type_ind)
                        if cur_com_type_ind-1 ==0
                            gui_data.cur_com.type = gui_data.community_types{numel(gui_data.community_types)};
                        else
                            gui_data.cur_com.type = gui_data.community_types{cur_com_type_ind-1};
                        end
                        graph_regen();
                        gui_data.init_X=X;
                        gui_data.plot_orig_graph_is_on = true;
                    end
                case 'add'
                    gui_data.cur_com.size = gui_data.cur_com.size+10;
                    graph_regen();
                    gui_data.init_X=X;
                    gui_data.plot_orig_graph_is_on = true;
                    gui_data.clustered_mat_inds = 1:size(X,1);
                case 'subtract'
                    if gui_data.cur_com.size(1)-10>0
                        gui_data.cur_com.size = gui_data.cur_com.size-10;
                        graph_regen();
                        gui_data.init_X=X;
                        gui_data.plot_orig_graph_is_on = true;
                        gui_data.clustered_mat_inds = 1:size(X,1);
                    end

                case 'c'
                    %                         cur_par = gui_data.cur_com;
                    %                         graph_regen();
                    % X= (gui_data.init_X(1:100,1:100));
                    X_mask = zeros(size(gui_data.init_X));
                    sz1 = gui_data.cur_com.size(1) + gui_data.cur_com.size(end)*gui_data.cur_com.is_bipartite;
                    sz2 = gui_data.cur_com.size(end) + gui_data.cur_com.size(1)*gui_data.cur_com.is_bipartite;
                    X_mask(gui_data.cur_com.loc(1)+gui_data.cur_com.offset+[1:sz1]-1,gui_data.cur_com.loc(2)+gui_data.cur_com.offset+[1:sz2]-1) = 1;
                    X=gui_data.init_X.*X_mask;
                    gui_data.plot_orig_graph_is_on = true;
                case 'r'
                    %                         cur_par = init_par;
                    %                         graph_regen();
                    X=gui_data.init_X;
                    gui_data.plot_orig_graph_is_on = true;
                case 'b'
                    gui_data.cur_com.is_bipartite = ~gui_data.cur_com.is_bipartite;
                    graph_regen();
                    gui_data.init_X=X;
                    gui_data.plot_orig_graph_is_on = true;
                    gui_data.clustered_mat_inds = 1:size(X,1);
                case 'd'
                    if ~isempty(gui_data.cur_com.prev_sibling) || ~isempty(gui_data.cur_com.Parent)
                        if ~isempty(gui_data.cur_com.prev_sibling)
                            gui_data.cur_com = gui_data.cur_com.prev_sibling;
                            gui_data.cur_com.next_sibling.delete;
                        elseif ~isempty(gui_data.cur_com.Parent)
                            tmp = gui_data.cur_com.Parent;
                            gui_data.cur_com.delete;
                            gui_data.cur_com = tmp;
                        end

                        graph_regen();
                        gui_data.init_X=X;
                        gui_data.plot_orig_graph_is_on = true;
                        gui_data.clustered_mat_inds = 1:size(X,1);
                    end
                case 'i' %inserts a 10x10 clique after the current community
                    gui_data.cur_com.insert('clique',10);
                    gui_data.cur_com = gui_data.cur_com.next_sibling;
                    graph_regen();
                    gui_data.init_X=X;
                    gui_data.plot_orig_graph_is_on = true;
                    gui_data.clustered_mat_inds = 1:size(X,1);
            end
        end
        switch evnt.Key
            case 'n'
                gui_data.noise_modification_mode_is_on = mod(gui_data.noise_modification_mode_is_on+1,3);
                gui_data.structure_modification_mode_is_on = false;
            case 's'
                gui_data.structure_modification_mode_is_on = ~gui_data.structure_modification_mode_is_on;
                gui_data.noise_modification_mode_is_on = false;
            case 'rightarrow'
                if ~isempty(gui_data.cur_com.next_sibling)
                    gui_data.cur_com = gui_data.cur_com.next_sibling;
                end
            case 'leftarrow'
                if ~isempty(gui_data.cur_com.prev_sibling)
                    gui_data.cur_com = gui_data.cur_com.prev_sibling;
                end
            case 'return'
                if ~isempty(gui_data.cur_com.Children)
                    gui_data.cur_com = gui_data.cur_com.Children(1);
                end
                %                         graph_mod_lvls(end+1) = 1;
            case 'backspace'
                if ~isempty(gui_data.cur_com.Parent)
                    gui_data.cur_com = gui_data.cur_com.Parent;
                end
            case 'tab'
                gui_data.clustered_mat_cur_cluster = 1;
                gui_data.clustered_mat_inds = gui_data.inds_all(1:gui_data.clustered_mat_cl_sizes(gui_data.clustered_mat_cur_cluster));
                gui_data.plot_clustered_mat_is_on = true;

                gui_data.cluster_selection_mode_is_on = true;
        end
    end
end
gui_data.plot_handles_is_on = true;
gui_plot();
