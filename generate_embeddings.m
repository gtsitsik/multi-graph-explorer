function [U,A,B,L,R,M,method_specific_vars,iters,time,obj_cur] = generate_embeddings(X,labels,alg,L_type_ind,R,M,thres,max_iters,alg_opts,print_type,mtimesx_exists)

if ~exist('labels','var')
    labels = [];
end

if ~exist('alg','var') || isempty(alg) 
    alg = 2;
end

if ~exist('L_type_ind','var')|| isempty(L_type_ind) 
    L_type_ind = 2;
end

if ~exist('R','var')
    R = [];
end

if ~exist('M','var')
    M = [];
end

if ~exist('thres','var')
    thres = [];
end

if ~exist('max_iters','var')
    max_iters = [];
end

if ~exist('alg_opts','var')
    alg_opts = [];
end

if ~exist('print_type','var')
    print_type = [];
end

if ~exist('mtimesx_exists','var')
    mtimesx_exists = exist('mtimesx')==3;
end

L=[];
if L_type_ind == 1
    if all(X == permute(X,[2,1,3]),'all')
        L = X;
    else % Naive Symmetrization
        L = (X+permute(X,[2,1,3]))/2;
    end
elseif max(L_type_ind == 2)
    assert(min(X(:))>=0, "Graph weights have to be non-negative")
    D = sum(X,2);
    D_pinv = zeros(size(D));
    D_pinv(abs(D)>eps) = 1./D(abs(D)>eps);
    if all(X == permute(X,[2,1,3]),'all')
        D_pinv_sqrt = sqrt(D_pinv);
        L = D_pinv_sqrt.*X.*permute(D_pinv_sqrt,[2 1 3]);
    else
        % TODO: Allow user to set teleportation probability
        h = 0.0; % Teleportation probability
        P = (1-h)*(D_pinv.*X)+h*(1-eye(size(X,1)))/(size(X,1)-1);
        S=[];
        for i = 1:size(X,3)
            %             [S(:,1,i),~] = eigs(sparse(P(:,:,i)'),1,1);
            %             S(:,1,i) = abs(S(:,1,i));
            % TODO: Give user option for dense or sparse eigendecomposition
            % Alternative solution using dense eigendecomposition
            [UU,LL]=eig(P(:,:,i)');
            [~,ii]=min(abs(diag(LL)-1));
            S(:,1,i) = abs(UU(:,ii));
        end
        S_sqrt = sqrt(S);
        S_pinv_sqrt = zeros(size(S_sqrt));
        S_pinv_sqrt(abs(S_sqrt)>eps) = 1./S_sqrt(abs(S_sqrt)>eps);
        L_tmp = S_sqrt.*P.*permute(S_pinv_sqrt,[2,1,3]);
        L = (L_tmp+permute(L_tmp,[2,1,3]))/2;
    end
end
switch alg
    case 2
        [U,A,B,W,iters,time,obj_cur] = comclus(L,R,M,thres,max_iters,alg_opts.beta,alg_opts.rho,alg_opts.thres_inner,print_type,mtimesx_exists);
        %             [~,max_St_inds]=max(B',[],2);
        %             B_max = [1:size(B,1)]'==max_St_inds';
        method_specific_vars.W=W;
    case {3,4}
        if alg_opts.structure=="true"
            % If R>=R_true and M>=M_true then a number of component groups 
            % equal to the true number of view structures is created, and each 
            % component group is assigned a number of components equal to the 
            % true number of communities in the corresponding view structure.
            % Any remaining components are assigned to component groups 
            % randomly.
            % If M<M_true, then component groups are removed randomly.
            % If R<R_true, then components are removed randomly.
            tmp = [];
            ii=1;
            for i= randperm(numel(labels))
                tmp = [tmp i*ones(1,numel(unique(labels{ii})))];
                ii = ii+1;
            end
            if R >= numel(tmp)
                tmp = [tmp(:,randperm(numel(tmp))) randi(M,1,R-numel(tmp))];
            else
                tmp = tmp(:,randperm(numel(tmp),R));            
            end
            B = double([1:M]'==tmp);
        elseif alg_opts.structure=="random"
            B = double([1:M]'==randi(M,1,R));
        else
            error("Incorrect structure specified for Symmetric Richcom")
        end
        switch alg
            case 3
                [U,A,B,iters,time,obj_cur] = symmetric_richcom(L,R,M,thres,max_iters,alg_opts.rho,B,print_type,mtimesx_exists);
                method_specific_vars=[];
            case 4
                [U,A,B,iters,time,obj_cur] = CMNC(L,R,M,thres,max_iters,alg_opts.delta,B,print_type,mtimesx_exists);
                method_specific_vars=[];
        end
end
