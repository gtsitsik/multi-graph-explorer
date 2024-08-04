function  params_str = params2str(params,type)
% PARAMS2STR Generate string describing experiment parameters.
%    params_str = PARAMS2STR(params) returns a string form of the
%    parameter tree params. 
%
%    params_str = PARAMS2STR(params,'filename') returns a shorter string.

if ~exist('type','var')
    params_str = helper(params,"params");
elseif type == "filename"
    % TODO: Needs to be generalized
    params_str = strjoin(fieldnames(params.embedding_method),'--')+"_"+...
        strjoin(string([params.R(1) params.R(end)-params.R(max([end-1,1])) params.R(end)]),'-')+"R_"+...
        strjoin(string(params.L_type_ind),'-')+"L_type_ind_"+...
        numel(params.sample)+"samples_"+...
        strjoin(fieldnames(params.clustering_method),'--')+"_"+...
        strjoin(params.column_normalization_type,'--')+"_"+...
        params.workers+"workers";
else
    error("Incorrect string type");
end
end

function  params_str = helper(params,prnt)
params_str="";
if isstruct(params)
    names = fieldnames(params);
    for i=1:numel(names)
        params_str = params_str + ...
            helper(params.(names{i}),prnt+"."+names{i});
    end
else
    try
        if ~ismissing(params)
            params_str = prnt+": ["+strjoin(string(params),' ')+"]";
        else
            params_str = prnt+": nan";
        end
    catch ex % TODO: Handle general case
        params_str = prnt+":";
    end
    params_str = params_str + "\n";
end
end