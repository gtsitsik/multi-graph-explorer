function [params_all,inds_all] = generate_combinations(params,print_is_on)

% If no input is given the default parameter set is used.
if ~exist('params','var') || isempty(params)
    params.embedding_method.ComClus.beta = [0.1,0.2];
    params.embedding_method.ComClus.rho = linspace(0,0.16,6);
    params.embedding_method.ComClus.thres_inner = [1e-6];
    params.embedding_method.Symmetric_Richcom.structure = "random";
    params.embedding_method.Symmetric_Richcom.rho = linspace(0,0.2,6);
    params.R = 6:8;
    params.clustering_method.maximum.nofield = nan;
    params.thres = [1e-6];
    params.max_iters = [1000];
    params.M = 3;
    params.L_type_ind = [1 2];
    params.sample = [1:1];
    params.workers = 16;
    params.clustering_method.kmeans.replicates = [5 8];
    params.clustering_method.kmeans.clusters_num = "3 2 2";
    params.clustering_method.kmeans.row_normalization_type =["none","unit"];
    params.column_normalization_type =["none","sqrtB","B"]  ;
end
if ~exist('print_is_on','var') || isempty(print_is_on)
    print_is_on = false;
end
prm.params = params;
[params_all,inds_all] = helper(prm,print_is_on );
fprintf("\n")
% If no input is given the default parameter combinations are also printed.
if nargin==0
    for i = 1:numel(params_all)
        params_str = params2str(params_all{i}.params);
        fprintf(params_str)
        disp("~~~~")
    end
end

end

function [combs, combs_inds] = helper(params,print_is_on)
root_name = string(fieldnames(params));
params = params.(root_name);
combs = cell(0);
combs_inds = cell(0);
combs_num_all=[];
msg = "";
% Retrieves all children of all children of the root node
f1 = fieldnames(params);
if ~isempty(f1)
    f2 =cell(1,numel(f1));
    for i=1:numel(f1)
        if isstruct(params.(f1{i}))
            f2{i} = fieldnames(params.(f1{i}));
        else
            f2{i} = params.(f1{i});
        end
        combs_num_all(i) = numel(f2{i});
    end
else
    % TODO: Handle case where root node has no children
    %     combs_num_all=1;
end
combs_num = prod(combs_num_all);

% For each combination of parameters it produces further combinations for 
% each parameter that is not a leaf node.
for i = 1:combs_num
    if root_name == "params" && print_is_on
        prev_msg_len = strlength(msg)-1;
        perc = round(i/combs_num*100);
        msg =  "["+repmat('.',1,round(3*perc/10))+repmat(' ',1,30-round(3*perc/10))+"]"+perc+"%%";
        fprintf(repmat('\b',1,prev_msg_len)+msg)
    end
    inds = cell(1,numel(combs_num_all));
    [inds{:}] = ind2sub(combs_num_all ,i);

    % when the last dimensions of combs_num_all are all 1, inds needs to be
    % manually augmented
    while numel(inds) < numel(combs_num_all)
        inds{end+1} = 1;
%         inds2(end+1)=1;
    end

    combs_num2_all = nan(1,numel(f1));
    ch_combs = cell(1,numel(f1));
    ch_combs_inds = cell(1,numel(f1));
    for j = 1:numel(f1)
        if isstruct(params.(f1{j}))
            prm = struct;
           
            %             if ~isempty(fieldnames(params.(f1{j}).(f2{inds{j}})))
            prm.(f2{j}{inds{j}}) = params.(f1{j}).(f2{j}{inds{j}});
            %             else
            %                 prm.(f2{inds{j}}).empty_field= nan;
            %             end
            [ch_combs{j},ch_combs_inds{j}] = helper(prm,print_is_on);
            combs_num2_all(j) = numel(ch_combs{j});
        else
            ch_combs{j} = params.(f1{j})(inds{j});
            ch_combs_inds{j} = {[]};
            combs_num2_all(j) = 1;
        end
    end
    combs_num2 = prod(combs_num2_all);
    
    for j = 1:combs_num2
        inds2 = cell(1,numel(combs_num2_all));
        [inds2{:}] = ind2sub(combs_num2_all,j);
        combs_tmp = struct;
        combs_inds_tmp = [];
        for k = 1:numel(f1)
            if isstruct(params.(f1{k}))
                s = string(fieldnames(ch_combs{k}{inds2{k}}));
                if ~isempty(s)
                    combs_tmp.(root_name).(f1{k}).(s) = ...
                        ch_combs{k}{inds2{k}}.(s);
                else
                    combs_tmp.(root_name).(f1{k}) = struct;
                end
            else
                combs_tmp.(root_name).(f1{k}) = ch_combs{k}(inds2{k});
            end
            combs_inds_tmp = [combs_inds_tmp inds{k} ch_combs_inds{k}{inds2{k}}];
        end
        combs{end+1} = combs_tmp;
        combs_inds{end+1} = combs_inds_tmp;
    end
end
end

