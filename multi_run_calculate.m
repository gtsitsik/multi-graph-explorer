function [clusterings,embeddings,duration,duration_sum,params,data_dir] = multi_run_calculate(params,...
    clusterings_save_is_on,embeddings_save_is_on,print_type)
% MULTI_RUN_CALCULATE Calculate multiple embedding and clustering tasks for
% all combinations of input parameters.
%    [clusterings,embeddings] = MULTI_RUN_CALCULATE(params) returns two
%    structs containing the calculated quantities related to the
%    clusterings and the embeddings for all combinations of parameters
%    specified in params.
%
%    [clusterings,embeddings] = MULTI_RUN_CALCULATE(params,
%    clusterings_save_is_on,embeddings_save_is_on) can be used to
%    automatically save the clusterings and embeddings structs to .mat
%    files. If any of them is a string of a path containing saved data,
%    then these data are loaded instead and no calculations take place, and
%    any input parameters corresponding to these data are ignored. By
%    default clusterings_save_is_on and embeddings_save_is_on are set to
%    false and therefore no .mat files are created.
%
%    [clusterings,embeddings] = MULTI_RUN_CALCULATE(params,
%    clusterings_save_is_on,embeddings_save_is_on,print_type) can be used
%    to print intermediate messages about the progress of calculations by
%    setting print_type to "basic". By default, print_type is set to "none"
%    which prints no messages.
mtimesx_exists = exist('mtimesx','file');

if ~exist('params','var') || isempty(params)
    params = struct;
end
params_msg = " is missing. The default value will be used.";
if ~isfield(params,'graph_tree') 
    warning("'graph_tree'"+params_msg)
    size_all= {
        {30,20,10},
        {50,10},
        {10,50}
        };
    noise_level =0.01;
    sparsity_level_all = linspace(0.85,0.99,5);
    for k = 1:numel(sparsity_level_all)
        params.graph_tree(k) = graph_tree_root;
%             tmp = rand(100,100,5);
%             tmp = (tmp+permute(tmp,[2,1,3]))/2;
%             params.graph_tree(k).Data = tmp;
        for i = 1:numel(size_all)
            params.graph_tree(k).Children(i).slices_num = 3;
            params.graph_tree(k).Children(i).noise_level = noise_level;
            params.graph_tree(k).Children(i).sparsity_level = sparsity_level_all(k);
            for j = 1:numel(size_all{i})
                params.graph_tree(k).Children(i).Children(j).type = 'clique';
                params.graph_tree(k).Children(i).Children(j).size = size_all{i}{j}*2;
            end
        end
%         params.graph_tree(k).labels = [];
        [~,params.graph_tree(k)]=create_graph(params.graph_tree(k));
    end
end

if ~isfield(params,'embedding_method')
    warning("'embedding_method'"+params_msg)
    %     params.embedding_method.Symmetric_Richcom.rho = [0:0.4:1];
    %     params.embedding_method.Symmetric_Richcom.structure = "random";
    params.embedding_method.ComClus.beta = [1];
    params.embedding_method.ComClus.rho = [0];
    params.embedding_method.ComClus.thres_inner = [1e-6];
end
if ~isfield(params,'thres')
    warning("'thres'"+params_msg)
    params.thres = [1e-6];
end
if ~isfield(params,'max_iters')
    warning("'max_iters'"+params_msg)
    params.max_iters = [1000];
end


% if ~isfield(params,'I') || ~isfield(params,'K')
%     X = create_graph(params.graph_tree);
%     if ~isfield(params,'I')
%         params.I = size(X,1);
%     end
%     if ~isfield(params,'K')
%         params.K = size(X,3);
%     end
% end

if ~isfield(params,'R')
    warning("'R'"+params_msg)
    params.R = [5];
end
if ~isfield(params,'M')
    warning("'M'"+params_msg)
    params.M = [3];
end
if ~isfield(params,'L_type_ind')
    warning("'L_type_ind'"+params_msg)
    params.L_type_ind = [2];
end
if ~isfield(params,'sample')
    warning("'sample'"+params_msg)
    params.sample = 1:5;
end
if ~isfield(params,'clustering_method')
    warning("'clustering_method'"+params_msg)
    params.clustering_method.kmeans.replicates = [1];
    params.clustering_method.kmeans.clusters_num = "3 2 2";
    params.clustering_method.kmeans.row_normalization_type = ["unit"];
    %     params.clustering_method.maximum.nofield = "nofield";
end
if ~isfield(params,'column_normalization_type')
    warning("'column_normalization_type'"+params_msg)
    params.column_normalization_type =["B"];
end
if ~isfield(params,'clustered_entity')
    warning("'clustered_entity'"+params_msg)
    params.clustered_entity =["nodes","views"];
end
if ~isfield(params,'clustering_measure')
    warning("'clustering_measure'"+params_msg)
    params.clustering_measure = ...
        ["NMI","ARI","AMI","silhouette_equal","silhouette_empirical"];
end
if ~exist('print_type','var') || isempty(print_type)
    print_type = "basic";
end

if print_type == "time"
    params.workers = 0;
elseif ~isfield(params,'workers')
    warning("'workers'"+params_msg)
    params.workers = feature('numcores');
end

if any(print_type == ["time","nothing","basic"])
    embedding_method_print_type = "nothing";
elseif print_type == "detailed"
    embedding_method_print_type = "basic";
end

my_dataqueue = [];
if all(print_type ~= ["nothing","time"])
    if mtimesx_exists
        disp("mtimesx will be used for faster computations")
    else
        disp("mtimesx was not found. Built-in MATLAB operations will be used instead")
    end
    disp("=============================================")
    my_dataqueue = parallel.pool.DataQueue;
    afterEach(my_dataqueue, @update_progress_stats);
end


if ~exist('clusterings_save_is_on','var') || isempty(clusterings_save_is_on) || print_type == "time"
    clusterings_save_is_on = false;
end

if ~exist('embeddings_save_is_on','var') || isempty(embeddings_save_is_on) || print_type == "time"
    embeddings_save_is_on = false;
end

if islogical(embeddings_save_is_on)
    if  embeddings_save_is_on || clusterings_save_is_on
        foldername = 'Experiments';
        if ~isdir(foldername)
            mkdir(foldername);
        end

        foldername2 = params2str(params,"filename");
        folder_filenames = string({dir(foldername).name});
        data_dir = "";
        if prod(folder_filenames ~= foldername2)
            data_dir = foldername+"/"+foldername2;
        else
            foldername2_copy_num = 1;
            while  sum(folder_filenames == foldername2+"_copy"+foldername2_copy_num)
                foldername2_copy_num = foldername2_copy_num+1;
            end
            data_dir = foldername+"/"+foldername2+"_copy"+foldername2_copy_num;
        end
    end
elseif any(embeddings_save_is_on == "Experiments/"+string({dir("Experiments").name}))
    % TODO: Now it is only allowed to be relative path. Make it more flexible.
    data_dir = embeddings_save_is_on;
    filenames = string({dir(data_dir).name})';
    filenames = filenames(contains(filenames,"embeddings"));
    embeddings = [];
    for i = 1:numel(filenames)
        disp("Loading embeddings data: "+round(i/numel(filenames)*100)+"%")
        l = load(data_dir+"/"+filenames(i));
        if isfield(l,'params_embeddings')
            params_embeddings = l.params_embeddings;
        else
            embeddings = [embeddings l.embeddings];
        end
    end

    % When existing embeddings are loaded, the input parameters related to
    % embedding generation are overwritten by the parameters of the loaded
    % embeddings.
    params_tmp = params;
    params = params_embeddings;
    params.clustering_method = params_tmp.clustering_method;
    params.column_normalization_type = params_tmp.column_normalization_type;
    params.clustering_measure = params_tmp.clustering_measure;
    params.clustered_entity = params_tmp.clustered_entity;
    params.workers = params_tmp.workers;

    embeddings_experiments_num = numel(embeddings);
    duration_embeddings = zeros(1,embeddings_experiments_num);
    total_real_time_start = tic;
else
    error("embeddings_save_is_on is not a valid directory")
end

% Create a parallel pool with the specified number of workers.
tmp = gcp('nocreate');
if params.workers>=1 && (isempty(tmp) || tmp.NumWorkers~=params.workers)
    delete(tmp)
    parpool(params.workers);
    disp("=============================================")
elseif params.workers==0
    delete(tmp)
end
% Suppresses warnings on client and all workers
if ~isempty(gcp('nocreate'))
    pctRunOnAll warning('off','stats:kmeans:MissingDataRemoved');
    pctRunOnAll warning('off','stats:kmeans:EmptyCluster');
    pctRunOnAll warning('off','stats:kmeans:EmptyClusterRep');
    pctRunOnAll warning('off','MATLAB:eigs:SigmaNearExactEig');
    pctRunOnAll warning on verbose
else
    warning('off','stats:kmeans:MissingDataRemoved');
    warning('off','stats:kmeans:EmptyCluster');
    warning('off','stats:kmeans:EmptyClusterRep');
    warning('off','MATLAB:eigs:SigmaNearExactEig');
    warning on verbose
end

% ~~~~~~~~~ Embeddings ~~~~~~~~~
embeddings_params = [];
if islogical(embeddings_save_is_on)
    params_embeddings = params;
    if isfield(params,'clustering_method')
        params_embeddings = rmfield(params_embeddings,'clustering_method');
    end
    if isfield(params,'column_normalization_type')
        params_embeddings = rmfield(params_embeddings,'column_normalization_type');
    end
    if isfield(params,'clustering_measure')
        params_embeddings = rmfield(params_embeddings,'clustering_measure');
    end
    if isfield(params,'clustered_entity')
        params_embeddings = rmfield(params_embeddings,'clustered_entity');
    end

    if all(print_type ~= ["nothing","time"])
        disp("Generating embedding parameters combinations")
        tic
    end
    embeddings_params = generate_combinations(params_embeddings,all(print_type ~= ["nothing","time"]));
    if all(print_type ~= ["nothing","time"])
        disp("Generation time: "+char(seconds(toc),'hh:mm:ss.SSSS'))
    end
    embeddings_experiments_num = numel(embeddings_params);

    % Shuffles experiments for more accurate predictions of remaining time
    % and better computational load distribution among parallel workers.
    embeddings_params = embeddings_params(randperm(embeddings_experiments_num));

    embeddings = struct(...
        'params_inds',cell(1,embeddings_experiments_num),...
        'data',cell(1,embeddings_experiments_num));
    if print_type ~= "time"
        % TODO: add message for completion percentage
        for i = 1:numel(embeddings)
            %     parfor(i = 1:numel(embeddings),params.workers)
            embeddings(i).params_inds = params2inds(embeddings_params{i}.params,params_embeddings);
        end
    end
    total_time = [];
    remaining_time = [];
    time_start = tic;
    cur_experiments_num = embeddings_experiments_num;
    num_iters_completed = 0;
    cur_data_type = "embeddings";
    prev_print_minute=0;
    X_time_all=cell(1,embeddings_experiments_num);
    graph_tree_time_all(1,1:embeddings_experiments_num) = graph_tree_root; 
    if print_type == "time"
        graph_tree_time_start = tic;
        for i = 1:embeddings_experiments_num
            cur_params = embeddings_params{i}.params;
            graph_tree = cur_params.graph_tree;
            [X_time_all{i},graph_tree_time_all(i)] = create_graph(graph_tree);
        end
        disp("Graph tree generation time: "+char(seconds(toc(graph_tree_time_start)),'hh:mm:ss'));
    end
    total_real_time_start = tic;

    duration_embeddings = nan(1,embeddings_experiments_num);
    %         for par_for_ind = 1:embeddings_experiments_num
    parfor(par_for_ind = 1:embeddings_experiments_num,params.workers)
        cur_params = embeddings_params{par_for_ind}.params;
        graph_tree = cur_params.graph_tree;
        R = cur_params.R;
        M = cur_params.M;
        L_type_ind = cur_params.L_type_ind;
        thres = cur_params.thres;
        max_iters = cur_params.max_iters;
        cur_embedding_method_name = string(fieldnames(cur_params.embedding_method));
        cur_alg_opts_names = fieldnames(cur_params.embedding_method.(cur_embedding_method_name));
        alg_opts = struct;
        for i=1:numel(cur_alg_opts_names)
            alg_opts.(cur_alg_opts_names{i})=cur_params.embedding_method.(cur_embedding_method_name).(cur_alg_opts_names{i});
        end

        if print_type ~= "time"
            [X,graph_tree] = create_graph(graph_tree);
        else
            X = X_time_all{par_for_ind};
            graph_tree = graph_tree_time_all(par_for_ind);
        end

        nodes_labels = {graph_tree.Children.labels};
        views_labels = graph_tree.labels;

        alg = [];
        switch cur_embedding_method_name
            case 'ComClus'
                alg = 2;
            case 'Symmetric_Richcom'
                alg = 3;
            case 'CMNC'
                alg = 4;
            otherwise
                error('Incorrect method name')
        end

        embeddings_real_time_start = tic;
        [U,A,B,L,R,M,embedding_method_vars,iters,fit_time,obj_cur] = generate_embeddings(X,nodes_labels,alg,L_type_ind,R,M,thres,max_iters,alg_opts,embedding_method_print_type,mtimesx_exists);
        duration_embeddings(par_for_ind) = toc(embeddings_real_time_start);

%         if any(isnan(U(:)))
%             disp("any(isnan(U(:)))")
%         end
%         if any(isnan(A(:)))
%             disp("any(isnan(A(:)))")
%         end
%         if any(isnan(B(:)))
%             disp("any(isnan(B(:)))")
%         end

        if alg == 2
            %         W_all{par_for_ind} = embedding_method_vars.W;
        end

        embeddings(par_for_ind).data.U_all=U;
        embeddings(par_for_ind).data.A_all=A;
        embeddings(par_for_ind).data.B_all=B;
        if print_type ~= "time"
            embeddings(par_for_ind).data.iters_all=iters;
            embeddings(par_for_ind).data.fit_times_all=fit_time;
            embeddings(par_for_ind).data.obj_cur_all=obj_cur;
            embeddings(par_for_ind).data.nodes_labels_all= nodes_labels;
            embeddings(par_for_ind).data.views_labels_all= views_labels;
            embeddings(par_for_ind).data.a_all = params2str(cur_params);
            if print_type ~= "nothing"
                send(my_dataqueue, NaN);
            end
        end
    end

    if embeddings_save_is_on
        save_results("embeddings");
    end
end
if all(print_type ~= ["nothing","time"])
    disp("=============================================")
end

% ~~~~~~~~~ Clusterings ~~~~~~~~~
params_cluster.clustering_method = params.clustering_method;
params_cluster.column_normalization_type = params.column_normalization_type;
params_cluster_2.clustering_measure = params.clustering_measure;
params_cluster_2.clustered_entity = params.clustered_entity;
if all(print_type ~= ["nothing","time"])
    disp("Generating clustering parameters combinations")
    tic
end
params_cluster_all = generate_combinations(params_cluster,all(print_type ~= ["nothing","time"]));
params_cluster_2_all = generate_combinations(params_cluster_2,all(print_type ~= ["nothing","time"]));
if all(print_type ~= ["nothing","time"])
    disp("Generation time: "+char(seconds(toc),'hh:mm:ss'))
end
clusterings_experiments_num = numel(params_cluster_all);
results_per_clustering_num = numel(params_cluster_2_all);
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
%==============================
% par_for_sizes = [embeddings_experiments_num clusterings_experiments_num];
% 
% clusterings = struct(...
%     'params',cell(1,prod(par_for_sizes)),...
%     'params_inds',cell(1,prod(par_for_sizes)),...
%     'data',cell(1,prod(par_for_sizes)));
% 
% 
% U_all=cellfun(@(x)x.U_all,{embeddings.data},'UniformOutput',false);
% A_all=cellfun(@(x)x.A_all,{embeddings.data},'UniformOutput',false);
% B_all=cellfun(@(x)x.B_all,{embeddings.data},'UniformOutput',false);
% if isempty(embeddings_params)
%     params_inds_all={embeddings.params_inds};
% end
% total_time = [];
% remaining_time = [];
% time_start = tic;
% cur_experiments_num = prod(par_for_sizes);
% num_iters_completed = 0;
% cur_data_type = "clusterings";
% prev_print_minute = 0;
% % for par_for_ind = 1:prod(par_for_sizes)
% parfor(par_for_ind = 1:prod(par_for_sizes),params.workers)
%     [cur_embeddings_ind,cur_cluster_ind] = ind2sub(par_for_sizes,par_for_ind);
%     U = U_all{cur_embeddings_ind};
%     A = A_all{cur_embeddings_ind};
%     B = B_all{cur_embeddings_ind};
%     
%     if ~isempty(embeddings_params)
%         cur_embeddings_params = embeddings_params{cur_embeddings_ind}.params;
%     else
%         cur_embeddings_params = inds2params(params_inds_all{cur_embeddings_ind},params_embeddings);
%     end
%     clusterings(par_for_ind).params = ...
%         mergestructs(...
%         cur_embeddings_params,...
%         params_cluster_all{cur_cluster_ind}.params...
%         );
%     % TODO: Add each clustering quality metric as parameter in the parameter tree
%     [tmp,~] = cluster_embeddings(U,A,B,clusterings(par_for_ind).params);
%     if print_type ~= "time"
%         clusterings(par_for_ind).data.nodes.cluster_qual = tmp.nodes.cluster_qual;
%         clusterings(par_for_ind).params_inds= params2inds(clusterings(par_for_ind).params,params);
%         if print_type ~= "nothing"
%             send(my_dataqueue, NaN);
%         end
%     end
% end
% if all(print_type ~= ["nothing","time"])
%     disp("=============================================")
% end
% disp("Average embedding and clustering time: "+char(seconds(toc(total_real_time_start)/numel(params.sample)),'hh:mm:ss.SSSS'));
% clusterings = rmfield(clusterings,'params');
%========================

if params.workers>0
    U_all=distributed(cellfun(@(x)x.U_all,{embeddings.data},'UniformOutput',false));
    A_all=distributed(cellfun(@(x)x.A_all,{embeddings.data},'UniformOutput',false));
    B_all=distributed(cellfun(@(x)x.B_all,{embeddings.data},'UniformOutput',false));
    duration_embeddings = distributed(duration_embeddings);
    if isempty(embeddings_params)
        params_inds_all = distributed({embeddings.params_inds});
    else
        embeddings_params = distributed(embeddings_params);
    end
else
    U_all=cellfun(@(x)x.U_all,{embeddings.data},'UniformOutput',false);
    A_all=cellfun(@(x)x.A_all,{embeddings.data},'UniformOutput',false);
    B_all=cellfun(@(x)x.B_all,{embeddings.data},'UniformOutput',false);
    if isempty(embeddings_params)
        params_inds_all={embeddings.params_inds};
    end
end

%FIXME: clear only when the function does not need to return the embeddings
% clear embeddings;

total_time = [];
remaining_time = [];
time_start = tic;
par_for_sizes = [embeddings_experiments_num clusterings_experiments_num];
cur_experiments_num = prod(par_for_sizes);
num_iters_completed = 0;
cur_data_type = "clusterings";
prev_print_minute = 0;
spmd(params.workers)
    if params.workers>0
        U_all_local = getLocalPart(U_all);
        A_all_local = getLocalPart(A_all);
        B_all_local = getLocalPart(B_all);
        duration_embeddings_local = getLocalPart(duration_embeddings);
        if isempty(embeddings_params)
            params_inds_all_local = getLocalPart(params_inds_all);
        else
            embeddings_params_local = getLocalPart(embeddings_params);
        end
        par_for_sizes_local = [numel(U_all_local),...
            clusterings_experiments_num,results_per_clustering_num];
    else
        U_all_local = U_all;
        A_all_local = A_all;
        B_all_local = B_all;
        duration_embeddings_local = duration_embeddings;
        if isempty(embeddings_params)
            params_inds_all_local = params_inds_all;
        else
            embeddings_params_local = embeddings_params;
        end
        par_for_sizes_local = [numel(U_all),clusterings_experiments_num,...
            numel(params_cluster_2_all)];

    end
    clusterings = struct(...
        'params_inds',cell(1,prod(par_for_sizes_local)),...
        'data',cell(1,prod(par_for_sizes_local)));

    duration_clusterings = nan(1,prod(par_for_sizes_local(1:2)));
    duration_sum = nan(1,prod(par_for_sizes_local(1:2)));
    % randperm is used to get more accurate estimations for the remaining
    % time of calculations
    for par_for_ind = randperm(prod(par_for_sizes_local(1:2)))
        [cur_embeddings_ind,cur_cluster_ind] = ...
            ind2sub(par_for_sizes_local(1:2),par_for_ind);
        U = U_all_local{cur_embeddings_ind};
        A = A_all_local{cur_embeddings_ind};
        B = B_all_local{cur_embeddings_ind};

        if ~isempty(embeddings_params)
            cur_embeddings_params = embeddings_params_local{cur_embeddings_ind}.params;
        else
            cur_embeddings_params = inds2params(params_inds_all_local{cur_embeddings_ind},params_embeddings);
        end

        cur_clusterings_params = ...
            mergestructs(...
            cur_embeddings_params,...
            params_cluster_all{cur_cluster_ind}.params...
            );
        cur_clusterings_params.workers = params.workers;
        clusterings_real_time_start = tic;
        [pred,~] = cluster_embeddings(U,A,B,cur_clusterings_params,print_type);
        duration_clusterings(par_for_ind) = toc(clusterings_real_time_start);
        duration_sum(par_for_ind) = duration_embeddings_local(cur_embeddings_ind)+duration_clusterings(par_for_ind);
        if print_type ~= "time"
            for i = 1:results_per_clustering_num
                ind = sub2ind(par_for_sizes_local,...
                    cur_embeddings_ind,cur_cluster_ind,i);
                tmp = params_cluster_2_all{i}.params;
                cur_clusterings_params_2 = mergestructs(...
                    cur_clusterings_params,tmp);
                if isfield(pred,tmp.clustered_entity)
                    tmp2 = pred.(tmp.clustered_entity).cluster_qual;
                else
                    tmp2 = [];
                end
                if isfield(tmp2,tmp.clustering_measure)
                    clusterings(ind).data = tmp2.(tmp.clustering_measure);
                else 
                    clusterings(ind).data = missing;
                end
                clusterings(ind).params_inds = params2inds(...
                    cur_clusterings_params_2,params);
            end
            if print_type ~= "nothing"
                send(my_dataqueue, NaN);
            end
        end
    end

end
duration_sum = [duration_sum{:}];
clusterings = [clusterings{:}];


if all(print_type ~= ["nothing","time"])
    disp("=============================================")
end
duration = toc(total_real_time_start);
disp("Average embedding and clustering time: "+char(seconds(duration/numel(clusterings)),'hh:mm:ss.SSSS'));

%===========================
if clusterings_save_is_on
    save_results("clusterings");
end
if all(print_type ~= ["nothing","time"])
    disp("=============================================")
end

% ~~~~~~~~~ Helper functions ~~~~~~~~~

    function update_progress_stats(data)
        num_iters_completed = num_iters_completed + 1;
        total_time(num_iters_completed) = toc(time_start);
        perc_completed = round(num_iters_completed/cur_experiments_num*100,0);
        num_iters_remaining = cur_experiments_num-num_iters_completed;
        iters_interv = 500;
        if num_iters_completed<=iters_interv
            remaining_time(num_iters_completed) = nan;
        else
            remaining_time(num_iters_completed) =  ((total_time(num_iters_completed)-total_time(num_iters_completed-iters_interv))/iters_interv)*num_iters_remaining;
        end
        cur_duration_sec = seconds(total_time(num_iters_completed));
        cur_print_minute = minutes(cur_duration_sec);

        % print every 15 seconds
        if floor(cur_print_minute*60/15)>prev_print_minute || num_iters_completed == cur_experiments_num
            disp(char(datetime('now'))+" | "+cur_data_type+": "+perc_completed+"% completed | time passed: "+ char(cur_duration_sec,'hh:mm:ss') +" | estimated remaining time:"+char(seconds(remaining_time(num_iters_completed)),'hh:mm'))
            prev_print_minute = cur_print_minute*60/15;
        end
    end

    function save_results(data_type)
        % Save results to a .mat file
        if ~isdir(data_dir)
            mkdir(data_dir);
        end
        cur_data =[];

        max_data_size = 2*1024^3-1;
        cur_start_ind = 1;
        cur_part = 0;
        filename =  data_type+"_part"+cur_part;

        switch data_type
            case "clusterings"
                save(char(data_dir+"/"+filename),'params','-v6');
                cur_data = clusterings;
            case "embeddings"
                save(char(data_dir+"/"+filename),'params_embeddings','-v6');
                cur_data = embeddings;
        end

        while cur_start_ind <= numel(cur_data)
            cur_part = cur_part+1;
            filename =  data_type+"_part"+cur_part;
            step = numel(cur_data)-cur_start_ind;
            cur_end_ind = cur_start_ind+step;
            while true
                data_type_tmp = cur_data(cur_start_ind:cur_end_ind);
                cur_data_size = whos('data_type_tmp').bytes;
                %                 num2str( [cur_part,cur_data_size,max_data_size, step])
                % TODO: This seems like an inefficient implementation of binary search where bytes are calculated about twice as many times than necessary.
                if step ==0
                    break
                elseif cur_data_size > max_data_size
                    cur_end_ind = cur_end_ind-step;
                    step = floor(step/2);
                    %                     data_type_tmp = cur_data(cur_start_ind:cur_end_ind);
                    %                     break;
                elseif cur_end_ind < numel(cur_data)
                    while cur_end_ind+step > numel(cur_data)
                        step = floor(step/2);
                    end
                    cur_end_ind = cur_end_ind+step;
                else
                    break;
                end
            end

            switch data_type
                case "clusterings"
                    clusterings = data_type_tmp;
                case "embeddings"
                    embeddings = data_type_tmp;
            end

            save(char(data_dir+"/"+filename),data_type,'-v6');
            cur_start_ind = cur_end_ind+1;
        end
        switch data_type
            case "clusterings"
                clusterings = cur_data;
            case "embeddings"
                embeddings = cur_data;
        end
        if print_type~="nothing"
            disp("Files saved");
        end
    end
end
