function [pred,U_nrm] = cluster_embeddings(U,A,B,params,print_type)
%CLUSTER_EMBEDDINGS Calculate view and node clusterings based on the
%embeddings of a multiview graph.
%
%    pred = CLUSTER_EMBEDDINGS(U,A,B) takes as input the view embeddings
%    defined as the rows of A, and the node embeddings of the i-th view
%    cluster defined as a function of U and the i-th row of B. It then
%    returns clustering information for the nodes and views in pred.nodes
%    and pred.views respectively. pred.nodes and pred.views are structured
%    very similarly and, therefore, below we only discuss the former:
%
%    pred.nodes.labels{i}        -  The node clustering for the i-th
%                                   calculated view cluster.
%    pred.nodes.cluster_qual.    -  Average Micro Silhouette Coefficient
%          silhouette_empirical     across all node clusterings.
%    pred.nodes.cluster_qual.    -  Average Macro Silhouette Coefficient
%              silhouette_equal     across all node clusterings.
%                                
%    
%    pred = CLUSTER_EMBEDDINGS(U,A,B,params) allows the user to specify the
%    parameters for the clustering as follows:
%
%    params.graph_tree.labels          -  Ground-truth view labels.
%    params.graph_tree.Children(i).    -  Ground-truth node labels for the
%                              labels     i-th view cluster.
%    params.clustering_method.         -  Value of k for applying k-means
%              kmeans.clusters_num(i)     on the i-th node clustering.
%    params.clustering_method.         -  The number of times k-means is
%                   kmeans.replicates     applied.
%    params.clustering_method.         -  Can have any value and is only
%                     maximum.nofield     for indicating that the maximum    
%                                         likelihood method should be used.
%    params.clustering_method.         -  Threshold (0-1) for the
%              large_inner_prod.thres     largest-inner-product method 
%    params.column_normalization_type  -  Normalization of the columns of
%                                         U.For the i-th node clustering
%                                         "B" multiplies the j-th column of
%                                         U with the j-th element of the
%                                         i-th row of B. "sqrtB" does the
%                                         same but modifies the elements of
%                                         B first as the square root of
%                                         their absolute value. "none"
%                                         performs no normalization.
%
%    If ground-truth view and node labels are specified, then extrinsic
%    measures for the clustering quality of nodes are also computed.
%    Specifically, we have
%
%    pred.nodes.cluster_qual.NMI(i)   -  Normalized Mutual Information for
%                                        the i-th ground-truth view
%                                        cluster.
%
%    pred.nodes.cluster_qual.ARI(i)   -  Adjusted Rand Index for the node
%                                        clustering corresponding to the
%                                        i-th ground-truth view cluster. 
%
%    pred.nodes.cluster_qual.AMI(i)   -  Adjusted Mutual Information for
%                                        the node clustering corresponding
%                                        to the i-th ground-truth view
%                                        cluster. 
%
%    NOTE: If the number of calculated view clusters is less than the
%    ground-truth ones, then the clustering performance for the additional
%    ground-truth view clusters is defined to be 0 for extrinsic measures.
%    If there are more calculated view clusters than ground-truth ones then
%    the additional ones are ignored in the calculations of extrinsic
%    measures.
% 
% 
%    pred = CLUSTER_EMBEDDINGS(...,print_type) shows basic information
%    about the progress if print_type is set to "basic" (default), while
%    more detailed information can be shown if it is set to "all". If it is
%    set to "time" then all messages are suppressed and no clustering
%    quality measure is computed.

%    -------------------------------------------------------------------------
%    set cl_num(i)=-1 for automatic estimation of communities in the i-th view
%    cluster
%    -------------------------------------------------------------------------
validate_input();

nodes_clustering_method = string(fieldnames(params.clustering_method));
if ~isfield(params,'graph_tree') || ...
        isempty(params.graph_tree.Children) ||  ( ...
        numel(unique(params.graph_tree.labels))==1 && ...
        all(cellfun(@(x)numel(unique(x))==1, ...
        {params.graph_tree.Children.labels})))
    labels_exist = false;
else
    labels_exist = true;
    views_labels = params.graph_tree.labels;
    nodes_labels = {params.graph_tree.Children.labels};
end
if nodes_clustering_method == "kmeans"
    tmp = params.clustering_method.kmeans.clusters_num;
    nodes_clusters_num = str2double(split(tmp)');
else
    nodes_clusters_num = nan(1,size(A,2));
end

[U_nrm,pred.nodes.labels]= deal(cell([1,size(A,2)]));

% Note that each column of A can be arbitrarily scaled by also scaling the
% corresponding row of B. This can lead to issues in ComClus and Symmetric
% Richcom where each row of A has more than one non-zero values.
% tmp = vecnorm(B,2,2);
% A = A.*tmp';
% B = B./tmp;

% rows of A should be clustered based on the value of maximum magnitude
cur_cluster_inds = my_cluster(A,"maximum");

A_cluster_inds_matrix = cur_cluster_inds==[1:max([cur_cluster_inds;size(A,2)])];

if labels_exist
    % ========= check - verify that prm is always a permutation matrix
    costt = normalize_fibers(A_cluster_inds_matrix,1)'*...
        normalize_fibers(views_labels'==unique(views_labels),1);
    [mtch,uR,~] = matchpairs(costt,0,'max');


    % unmatched clusters are reordered so that clusters with lower median
    % slice_id have lower pred.nodes.labels.
    if ~isempty(uR)
        tmp = A_cluster_inds_matrix(:,uR).*[1:size(A_cluster_inds_matrix,1)]';
        tmp(tmp==0)=nan;
        [~,tmp] = sort(median(tmp,1,'omitnan'));
        tmp2 = setdiff([1:size(costt,1)]',mtch(:,2));
        mtch = [mtch; uR(tmp) tmp2(1:numel(uR))];
    end


    prm = zeros(size(costt,1),max(size(costt)));
    prm(sub2ind(size(prm),mtch(:,1),mtch(:,2))) = 1;
    prm_small = prm(1:size(A,2),:);

    A = A*prm_small;
    A_cluster_inds_matrix = A_cluster_inds_matrix*prm;
    B = prm_small'*B;
end

B_useful_inds = abs(B)>=eps;
for i = 1:size(A,2)
    tmp = B_useful_inds(i,:);
    switch params.column_normalization_type
        case "none"
            U_nrm_tmp = U(:,tmp);
        case "sqrtB"
            U_nrm_tmp = U(:,tmp).*sqrt(abs(B(i,tmp)));
        case "B"
            U_nrm_tmp = U(:,tmp).*B(i,tmp);
        otherwise
            tmp = params.column_normalization_type;
            error("'"+tmp+"' is not a valid column normalization type")
    end

    switch nodes_clustering_method
        case "kmeans"
            switch params.clustering_method.kmeans.row_normalization_type
                case "none"
                    U_nrm{i} = U_nrm_tmp;
                case "unit"
                    %nodes with small norm are not clustered
                    U_nrm{i} = normalize_fibers(U_nrm_tmp,2);
                otherwise
                    tmp = params.clustering_method.kmeans;
                    tmp = tmp.row_normalization_type;
                    error("'"+tmp+"' is not a valid row normalization type")
            end

            % automatic number of clusters applies only for k-means
            if i>numel(nodes_clusters_num)
                nodes_clusters_num(i) = 0;
            elseif nodes_clusters_num(i) == -1
                nodes_clusters_num(i) = size(U_nrm{i},2);
            end
        case "maximum"
            U_nrm{i} = U_nrm_tmp;
            nodes_clusters_num(i) = size(U_nrm{i},2);
        case "large_inner_prod"
            U_nrm{i} = U_nrm_tmp;
            nodes_clusters_num(i) = nan;
        otherwise
            tmp = nodes_clustering_method;
            error("'"+tmp+"' is not a valid clustering method")

    end
end


% A_cluster_inds_matrix(:,sum(A_cluster_inds_matrix,1)==0)=[];
pred.views.labels = ...
    sum((A_cluster_inds_matrix).*[1:size(A_cluster_inds_matrix,2)],2)';
uniq_pred_view_labels = unique(pred.views.labels);

if print_type ~= "time"
    [silhouette_equal_all,silhouette_empirical_all,...
        NMI_all,ARI_all,AMI_all]=deal([]);
    if labels_exist
        qual_metrics = py.importlib.import_module('sklearn.metrics');
        py_nmi = @(x,y) qual_metrics.normalized_mutual_info_score(x,y);
        py_ari = @(x,y) qual_metrics.adjusted_rand_score(x,y);
        py_ami = @(x,y) qual_metrics.adjusted_mutual_info_score(x,y);
        [NMI_all,AMI_all,ARI_all] = deal(zeros(1,numel(nodes_labels)));
    end
    [silhouette_empirical_all,silhouette_equal_all] = ...
        deal(zeros(1,numel(uniq_pred_view_labels)));
end

for j = 1:numel(uniq_pred_view_labels)
    i = uniq_pred_view_labels(j);
    if nodes_clusters_num(i)==0 
        continue
    end

    nodes_clustering_opts = params.clustering_method.(nodes_clustering_method);
    nodes_clustering_opts.clusters_num = nodes_clusters_num(i);

    cur_cluster_inds = my_cluster(U_nrm{i},...
        nodes_clustering_method,nodes_clustering_opts);
    if nodes_clustering_method == "large_inner_prod"
        nodes_clusters_num(i) = numel(unique(cur_cluster_inds));
    end

    pred.nodes.labels{j} = cur_cluster_inds;
    if print_type ~= "time"
        % TODO: Allow for non-numeric nodes_labels
        if labels_exist && i<=numel(nodes_labels)
            NMI_all(i) = py_nmi(cur_cluster_inds,nodes_labels{i});
            ARI_all(i) = py_ari(cur_cluster_inds,nodes_labels{i});
            AMI_all(i) = py_ami(cur_cluster_inds,nodes_labels{i});
        end

        if  size(U_nrm{i},2)~=0 %&& ~all(abs(A(:,i))<eps)
            % TODO: Check why silhouette produces many NaN's
            tmp = silhouette(U_nrm{i},cur_cluster_inds)';
            if print_type == "all"
                if size(U_nrm{i},2)~=0 && ( all(abs(A(:,i))<eps) ||  all(abs(B(i,:))<eps) )
                    disp("size(U_nrm{i},2)~=0 && ( all(abs(A(:,i))<eps) ||  all(abs(B(i,:))<eps) )");
                end
                if any(isnan(tmp)) && ~all(isnan(tmp))
                    disp("any(isnan(tmp)) && ~all(isnan(tmp)) | "+numel(unique(cur_cluster_inds)))
                end
            end
            silhouette_empirical_all(j) = mean(tmp,'omitnan');

            tmp2 = [];
            for k = unique(cur_cluster_inds')
                tmp2(end+1) = mean(tmp(cur_cluster_inds==k));
            end
            silhouette_equal_all(j) = mean(tmp2,'omitnan');
        end
    end
end

if print_type ~= "time"
    pred.nodes.cluster_qual.silhouette_empirical =mean(silhouette_empirical_all,'omitnan');
    pred.nodes.cluster_qual.silhouette_equal = mean(silhouette_equal_all,'omitnan');
    view_l_pred = pred.views.labels;
    tmp = silhouette(A,view_l_pred)';
    pred.views.cluster_qual.silhouette_empirical = mean(tmp,'omitnan');
    tmp2 = [];
    for k = unique(view_l_pred)
        tmp2(end+1) = mean(tmp(view_l_pred==k));
    end
    pred.views.cluster_qual.silhouette_equal = mean(tmp2,'omitnan');
    if labels_exist
        pred.nodes.cluster_qual.NMI = NMI_all;
        pred.nodes.cluster_qual.ARI = ARI_all;
        pred.nodes.cluster_qual.AMI = AMI_all;
        %TODO: fix bug where NMI ARI and AMI throw error when only one view
        %exists
        if numel(views_labels)>1
            pred.views.cluster_qual.NMI = py_nmi(views_labels,view_l_pred);
            pred.views.cluster_qual.ARI = py_ari(views_labels,view_l_pred);
            pred.views.cluster_qual.AMI = py_ami(views_labels,view_l_pred);
        end
    end
end
% pred.nodes.clusters_num = nodes_clusters_num*prm_small';
pred.nodes.clusters_num = nodes_clusters_num;


    function validate_input()
        if ~exist('U','var') || isempty(U) || ...
                ~exist('A','var') || isempty(A) || ...
                ~exist('B','var') || isempty(B)
            err_ID = 'cluster_embeddings:not_enough_arguments';
            msg = 'Not enough input arguments';
            error(err_ID,msg);
        end
        if ~exist('print_type','var') || isempty(print_type)
            print_type = "basic";
        end
        default_params.clustering_method.large_inner_prod.thres = 0.8;
        default_params.column_normalization_type = "sqrtB";
        if ~exist('params','var') || isempty(params)
            params = default_params;
        else
            if ~isfield(params,'clustering_method')
                params.clustering_method = default_params.clustering_method;
            elseif numel(fieldnames(params.clustering_method))>1
                err_ID = 'cluster_embeddings:too_many_clustering_methods';
                msg = 'More than one clustering methods were specified.';
                error(err_ID,msg);
            elseif all(string(fieldnames(params.clustering_method)) ~= ...
                    ["kmeans","maximum","large_inner_prod"])
                err_ID = 'cluster_embeddings:invalid_clustering_method';
                msg = '"%s" is not a valid clustering method.';
                error(err_ID,msg,char(fieldnames(params.clustering_method)));
            end
            if ~isfield(params,'column_normalization_type')
                params.column_normalization_type = ...
                    default_params.column_normalization_type;
            elseif all(params.column_normalization_type ~= ["none","sqrtB","B"]);
                err_ID = 'cluster_embeddings:invalid_column_normalization_type';
                msg = '"%s" is not a valid column normalization type';
                error(err_ID,msg,char(params.column_normalization_type));
            end
        end
    end
end

function  cluster_inds = my_cluster(P,method_name,opts)
if isempty(P)
    cluster_inds = nan(size(P,1),1);
else
    switch method_name
        case "kmeans"
            try
                cluster_inds = kmeans(P,...
                    opts.clusters_num,...
                    'Replicates',opts.replicates,...
                    'EmptyAction','drop');
            catch exception
                % TODO: Check if this is still valid.
                % TODO: Handle warning about rows with missing data
                % for when there are less non-zero datapoints than clusters
                cluster_inds = nan(size(P,1),1);
%		disp(exception)
            end
        case "maximum" % maximum absolute value
            [cur_P_max_values,cluster_inds] = max(abs(P),[],2);
            cluster_inds(isnan(cur_P_max_values)) = nan;
        case "large_inner_prod"
            P = normalize_fibers(P,2);
            thres = opts.thres;
            structure = P*P';
            min_S = min(structure(:));
            max_S = max(structure(:));
            structure = (structure-min_S)/(max_S-min_S);
            cluster_inds = conncomp(graph(structure>thres))';
        otherwise
            error("Incorrect clustering type")
    end
end
% All samples that include NaN's are assigned to their own cluster which 
% contains no other samples.
tmp = isnan(cluster_inds);
cluster_inds(tmp) = max([cluster_inds; size(P,2)]) + [1:sum(tmp)];
end
