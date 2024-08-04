function prms = inds2params(prms_ind,params)
[prms,prms_ind] = helper(prms_ind,params);
if ~isempty(prms_ind)
    error('inds2params:tooManyInds','More than necessary number of indices')
end

end


function [prms,prms_ind] = helper(prms_ind,params)
f = fieldnames(params);
for i = 1:numel(f)
    if isempty(prms_ind)
        error('inds2params:insufficientInds','Insufficient number of indices')
    end
    if isstruct(params.(f{i}))
        f2 = fieldnames(params.(f{i}));
        if prms_ind(1) > numel(f2)
            error('inds2params:outOfRange','One or more indices are larger than the total number of the corresponding values')
        end
        tmp = f2{prms_ind(1)};
        prms_ind(1) = [];
        [prms.(f{i}).(tmp),prms_ind] = helper(prms_ind,params.(f{i}).(tmp));
    else
        if prms_ind(1) > numel(params.(f{i}))
            error('inds2params:outOfRange','One or more indices are larger than the total number of the corresponding values')
        end
        prms.(f{i}) = params.(f{i})(prms_ind(1));
        prms_ind(1) = [];
    end
end
end