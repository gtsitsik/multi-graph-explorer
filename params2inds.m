function [prms_ind,prms_names] = params2inds(prms,params)
[prms_ind,prms_names] = helper(prms,params,"");
end

% TODO: Remove first space

function [prms_ind,prms_names] = helper(prms,params,prefix)
if prefix~=""
    prefix = prefix+".";
end
f = fieldnames(params);
prms_ind = [];
prms_names = strings(0);
for i = 1:numel(f)
    if ~isfield(prms,f{i})
        error('params2inds:missing_field',"The given combination of parameters has a missing field")
    end
    if isstruct(params.(f{i}))
        f2 = fieldnames(params.(f{i}));
        is_valid = false;
        for j = 1:numel(f2)
            if  f2{j}==string(fieldnames(prms.(f{i})))
                [a,b]=helper(prms.(f{i}).(f2{j}),params.(f{i}).(f2{j}),prefix+string(f{i})+"."+string(f2{j}));
                prms_ind = [prms_ind, j, a];
                prms_names = [prms_names , prefix+string(f{i}), b];
                is_valid = true;
                break;
            end
        end
        if ~is_valid
            error('params2inds:invalid_value',"The given combination of parameters does not belong in the given range")
        end
    else
        f2 = params.(f{i});
        f3 = prms.(f{i});
        is_valid = false;
        for j = 1:numel(f2)
            % FIXME: Using simple handles comparison gives error probably
            % with graph_tree handles which seems to happen only on
            % parallel workers. Probably since handle objects are not the
            % same objects anymore when transfered to workers.

            % Comparing the handle objects with == can be much faster than
            % using isequal because it just compares their handles.
            % However, if two handle objects have identical properties but
            % different handles they will be considered to be different.
%             if isa(f2(j),'handle') && isa(f3,'handle') 
%                 tmp = f2(j)==f3;
%             else
                  tmp = isequal(f2(j),f3);
%             end
            if tmp
                prms_ind = [prms_ind j];
                prms_names = [prms_names , prefix+string(f{i})];
                is_valid = true;
                break;
            end
        end
        if ~is_valid
            error('params2inds:invalid_value',"The given combination of parameters does not belong in the given range")
        end
    end
end

end
