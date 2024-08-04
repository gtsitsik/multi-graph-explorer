function X_nrm = normalize_fibers(X,dim,small_norm_handling)
% NORMALIZE_FIBERS Normalize fibers of a tensor to be of unit 2-norm. 
%
%    X_nrm = NORMALIZE_FIBERS(X) returns a tensor X_nrm which is X with its
%    fibers across its first dimension normalized to be of unit 2-norm.
%
%    X_nrm = NORMALIZE_FIBERS(X,dim) normalizes the fibers across the
%    selected dimension dim.
%
%    X_nrm = NORMALIZE_FIBERS(X,dim,small_norm_handling) also specifies how
%    to handle fibers whose norm is very small. If small_norm_handling is:
%           "none" - no special handling
%            "NaN" - fibers with norm less than eps have 
%                    all their elements set to NaN

if nargin==1
    dim=1;
end
if nargin==2
    small_norm_handling = "none";
end
norm_X = vecnorm(double(X),2,dim);
if small_norm_handling == "none"
    norm_X(norm_X <= eps) = 1;
    X_nrm = X./norm_X;
elseif small_norm_handling == "NaN"
    norm_X(norm_X <= eps) = NaN;
    X_nrm = X./norm_X;
else
    error(small_norm_handling + " is not a valid value for small_norm_handling")
end

end
