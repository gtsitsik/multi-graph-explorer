function [X,par] = create_graph(par)
% TODO: Add documentation


% TODO: Allow the creation of real-world graphs with just the data and
% labels as direct inputs.
% TODO: Add global noise and sparsity to the graph_tree_root
% TODO: Add input validation checks for artificial graphs where the size
% and the type of each community need to be specified along with the number
% of slices.
% each community
if nargin == 0 % default graph. runs only the first time
    sizes= [20 20 20;30 10 0]*2;
    %     sizes= [30 20 10;50 10 0];
    % node_weights = 1;
    par = graph_tree_root;
    par.Children(1).is_symmetric = true;
    par.Children(1).slices_num = 3;
    par.Children(1).noise_level = 0.02;
    par.Children(1).sparsity_level = 0.85;
    par.Children(1).Children(1).size = sizes(1,1);
    par.Children(1).Children(1).type = 'clique';
    par.Children(1).Children(2).size = sizes(1,2);
    par.Children(1).Children(2).type = 'clique';
    par.Children(1).Children(3).size = sizes(1,3);
    par.Children(1).Children(3).type = 'clique';
    
    par.Children(2).is_symmetric = true;
    par.Children(2).slices_num = 1;
    par.Children(2).noise_level = 0.02;
    par.Children(2).sparsity_level = 0.85;
    par.Children(2).Children(1).size = sizes(2,1);
    par.Children(2).Children(1).type = 'clique';
    par.Children(2).Children(2).size = sizes(2,2);
    par.Children(2).Children(2).type = 'clique';
    par.Children(2).Children(3).size = sizes(2,3);
    par.Children(2).Children(3).type = 'clique';
    
    
end
if class(par)=="graph_tree_root" % runs only the first time
    data_is_given = ~isempty(par.Data);
    
    X=[];
    if data_is_given
        % If either view or node labels are not given, all views are added
        % to the same cluster, similarly all nodes are added to the same
        % cluster.
        if isempty(par.labels) || isempty(par.Children)
            par.labels = ones(1,size(par.Data,3));
            par.Children(1).labels = ones(1,size(par.Data,1));
        end
        unique_view_labels = unique(par.labels);
    else
        par.labels = []; % initialization

        % Augments each view structure with a dummy blank node cluster that
        % contains all the nodes found exclusively in the other view
        % structures.
        slice_size_all = zeros(1,numel(par.Children));
        dummy_all = {};
        for i = 1:numel(par.Children)
            dummy_all{i} = [];
            for ch = par.Children(i).Children
                if strcmp(ch.type,'dummy_blank')
                    dummy_all{i} = ch;
                else
                    slice_size_all(i) = slice_size_all(i) + ch.size;
                end
            end
            if isempty(dummy_all{i})
                dummy_all{i} = graph_tree_node;
                dummy_all{i}.type = 'dummy_blank';
                par.Children(i).Children(end+1) = dummy_all{i};
            end
        end
        max_size = max(slice_size_all);
        for i = 1:numel(dummy_all)
            dummy_all{i}.size = max_size - slice_size_all(i);
        end
    end

    for i=1:numel(par.Children)
        if data_is_given
            par.Children(i).type = 'given';
            par.Children(i).slices_num = numel(find(par.labels==unique_view_labels(i)));
            par.Children(i).size = size(par.Data,[1]);
        end
        cur_X = create_graph(par.Children(i));
        if ~data_is_given
            X(1:size(cur_X,1),1:size(cur_X,2),end+[1:size(cur_X,3)])=cur_X;
            par.labels = [par.labels i*ones(1,size(cur_X,3))];
        end
    end

    if ~data_is_given
        X(:,:,1)=[]; % first slice is always all-zeros for generated data
    else
        X = par.Data;
    end
else
    
    if isempty(par.Parent)
        par.loc = [1 1];
    end
    
    X=[];
    N_mask=[];
    
    switch par.type
        case 'given'
            % TODO: add sparsity, noise and modification for non-generated data
            
            %             unique_node_labels = unique(par.labels);
            %             for i = 1:numel(unique_node_labels)
            %                 par.Children(i).Parent = par;
            %                 par.Children(i).type = 'given (leaf)';
            %                 par.Children(i).slices_num = par.slices_num;
            %
            %
            %                 if i == 1
            %                     par.Children(i).loc = par.loc;
            %                     if numel(unique_node_labels)>1
            %                         par.Children(i+1) = graph_tree_node;
            %                         par.Children(i).next_sibling = par.Children(i+1);
            %                     end
            %                 else
            %                     par.Children(i).loc = par.Children(i-1).loc + par.Children(i-1).size(1)+par.Children(i-1).size(end)*par.Children(i-1).is_bipartite + par.Children(i).offset;
            %                     if i<numel(unique_node_labels)
            %                         par.Children(i+1) = graph_tree_node;
            %                         par.Children(i).next_sibling = par.Children(i+1);
            %                     end
            %                     par.Children(i).prev_sibling = par.Children(i-1);
            %                 end
            %
            %                 par.Children(i).labels = unique_node_labels(i)*ones(1,numel(find(unique_node_labels(i)==par.labels)));
            %
            %             end
            return
        case ''
            par.labels = []; % initialization
            for i = 1:numel(par.Children)
                par.Children(i).Parent = par;
                %             par.Children(i).type2 = par.type2;
                par.Children(i).slices_num = par.slices_num;
                
                if par.is_symmetric
                    par.Children(i).is_symmetric = true;
                end
                
                if i == 1
                    par.Children(i).loc = par.loc;
                    if numel(par.Children)>1
                        par.Children(i).next_sibling = par.Children(i+1);
                    end
                else
                    par.Children(i).loc = par.Children(i-1).loc + par.Children(i-1).size(1)+par.Children(i-1).size(end)*par.Children(i-1).is_bipartite + par.Children(i).offset;
                    if i<numel(par.Children)
                        par.Children(i).next_sibling = par.Children(i+1);
                    end
                    par.Children(i).prev_sibling = par.Children(i-1);
                end
                
                
                X_ch = create_graph(par.Children(i));
                if ~isempty(par.labels)
                    par.Children(i).labels = max(par.labels)+par.Children(i).labels;
                end
                par.labels = [par.labels par.Children(i).labels];
                %             if ~isempty(par.Children(i).type)
                %                 switch par.Children(i).slices_unique_entries_type
                %                     case ''
                % %                         if par.Children(i).is_symmetric
                % %                             inds = find(triu(X_ch_init,1)>0)';
                % %                         else
                % %                             inds = find(X_ch_init>0)';
                % %                         end
                %                         X_ch=[];
                %                         for j = 1:par.Children(i).slices_num'
                % %                             cur_inds = inds(randperm(numel(inds),ceil(rand*numel(inds))));
                %                             X_ch_tmp = X_ch_init(:,:,j);
                %
                % %                             X_ch_tmp(inds)=~X_ch_tmp(inds);
                %                             X_ch(:,:,j)=X_ch_tmp;
                %                         end
                %                     case 'equal'
                %                         if par.Children(i).is_symmetric
                %                             inds = find(triu(X_ch_init,1)>0)';
                %                         else
                %                             inds = find(X_ch_init>0)';
                %                         end
                %                         inds = inds(randperm(numel(inds),numel(inds)));
                %                         common_inds = inds(1:(ceil(numel(inds)*par.Children(i).slices_common_entries_perc)));
                %                         unique_inds = inds(numel(common_inds)+1:end) ;
                %                         tmp1 = numel(unique_inds);
                %                         tmp2 = par.Children(i).slices_num;
                %
                %                         if mod(tmp1,tmp2)>0
                %                             unique_inds = [unique_inds nan(1,tmp2-mod(tmp1,tmp2))];
                %                             unique_inds = unique_inds(randperm(numel(unique_inds),numel(unique_inds)));
                %                         end
                %                         unique_inds = reshape(unique_inds,[],par.Children(i).slices_num)';
                %                         X_ch=[];
                %                         for j = 1:par.Children(i).slices_num'
                %                             cur_inds = [common_inds unique_inds(j,:)];
                %                             cur_inds(isnan(cur_inds))=[];
                %                             X_ch_tmp = zeros(size(X_ch_init));
                %                             X_ch_tmp(cur_inds) = X_ch_init(cur_inds);
                %
                %                             %pads array in case various structures do not
                %                             %have equal size
                %                             X_ch(:,:,j) = padarray(X_ch_tmp,max([size(X_ch,1)-size(X_ch_tmp,1),0]),max([size(X_ch,2)-size(X_ch_tmp,2),0]),'post');
                % %                          X_ch(:,:,j) = X_ch_tmp;
                %
                %                         end
                %                     case 'random'
                %                         if par.Children(i).is_symmetric
                %                             X_ch_init_tmp = [];
                %                             for i = 1:size(X_ch_init,3)
                %                                 X_ch_init_tmp(:,:,i) = triu(X_ch_init(:,:,i),1);
                %                             end
                %                             inds = find(X_ch_init_tmp~=0)';
                %                         else
                %                             inds = find(X_ch_init>0)';
                %                         end
                %                         X_ch=[];
                %                         for j = 1:par.Children(i).slices_num'
                %                             cur_inds = inds(randperm(numel(inds),ceil(rand*numel(inds))));
                %                             X_ch_tmp = X_ch_init;
                %
                %                             X_ch_tmp(cur_inds)=~X_ch_tmp(cur_inds);
                %                             X_ch(:,:,j)=X_ch_tmp;
                %                         end
                %                 end
                %             else
                %                 X_ch = X_ch_init;
                %             end
                
                if ~par.Children(i).is_bipartite
                    I1 = [1:size(X_ch,1)]+size(X,1)+par.Children(i).offset;
                    I2 = [1:size(X_ch,2)]+size(X,2)+par.Children(i).offset;
                    X = padarray(X,[max(I1)-size(X,1) max(I2)-size(X,2)  ],'post');
                    X(I1,I2,1:size(X_ch,3)) = X_ch;
                    N_mask = padarray(N_mask,[max(I1)-size(N_mask,1) max(I2)-size(N_mask,2)  ],1,'post');
                    N_mask(I1,I2)= par.Children(i).has_global_noise+par.Children(i).has_global_sparsity*1j;
                else
                    I1 = [1:size(X_ch,1)]+size(X,1)+par.Children(i).offset;
                    I2 = [1:size(X_ch,2)+size(X_ch,1)]+size(X,2)+par.Children(i).offset;
                    X = padarray(X,[max(I2)-size(X,1) max(I2)-size(X,2)  ],'post');
                    X(I1,I2)= double(X(I1,I2) + [zeros(size(X_ch,1)) X_ch]);
                    N_mask = padarray(N_mask,[max(I2)-size(N_mask,1) max(I2)-size(N_mask,2)  ],1,'post');
                    N_mask(I1,I2)= par.Children(i).has_global_noise+par.Children(i).has_global_sparsity*1j;
                end
                
            end
        otherwise
            switch par.type
                case 'clique'
                    X_tmp = ones(par.size)-eye(par.size);
                case 'star'
                    X_tmp = zeros(par.size);
                    X_tmp([1:min(par.size)-1],[2:min(par.size)]) = diag(ones(1,min(par.size)-1));
                    X_tmp(1,3:end) = ones(1, par.size(end)-2);
                case 'triangular'
                    X_tmp = flip(tril(ones(par.size)));
                case 'diagonal'
                    X_tmp = ones(par.size);
                case {'blank','dummy_blank'}
                    X_tmp = zeros(par.size);
                otherwise
                    disp("Invalid community type")
            end
            X = repmat(X_tmp,1,1,par.slices_num);
    end
    % if ~isempty(par.type) %&& par.is_symmetric
    %     X = triu(X,1);
    % end
    if isempty(N_mask)
        N_mask = ~eye(size(X,[1,2]))*(1+1j);
    else
        N_mask = N_mask-diag(diag(N_mask));
    end
    
    % adding sparsity
    if par.is_symmetric
        X_tmp = [];
        for i = 1:size(X,3)
            X_tmp(:,:,i) = triu(X(:,:,i),1);
        end
        sparsity_inds_all = reshape(find(X_tmp.*imag(N_mask)~=0),1,[]);
        sparse_count = numel(sparsity_inds_all);
        sparsity_inds_inds = randperm(sparse_count,floor(sparse_count*par.sparsity_level));
        sparsity_inds_tmp = sparsity_inds_all(sparsity_inds_inds);
        [I,J,K] = ind2sub(size(X),sparsity_inds_tmp);
        sparsity_inds = [sparsity_inds_tmp,sub2ind(size(X),J,I,K)];
    else
        sparsity_inds_all = reshape(find(X.*imag(N_mask)~=0),1,[]);
        sparse_count = numel(sparsity_inds_all);
        sparsity_inds_inds = randperm(sparse_count,floor(sparse_count*par.sparsity_level));
        sparsity_inds = sparsity_inds_all(sparsity_inds_inds);
    end
    X(sparsity_inds)=0;

    
    % adding noise
    if par.is_symmetric
        N_mask_tmp = [];
        for i = 1:size(X,3)
            N_mask_tmp(:,:,i) = triu(N_mask,1);
        end
        noise_inds_all = reshape(find(real(N_mask_tmp)~=0),1,[]);
        noise_count = numel(noise_inds_all);
        noise_inds_inds = randperm(noise_count,floor(noise_count*par.noise_level));
        noise_inds_tmp = noise_inds_all(noise_inds_inds);
        [I,J,K] = ind2sub(size(X),noise_inds_tmp);
        noise_inds = [noise_inds_tmp,sub2ind(size(X),J,I,K)];
    else
        noise_inds_all = reshape(find(repmat(real(N_mask),1,1,size(X,3))~=0),1,[]);
        noise_count = numel(noise_inds_all);
        noise_inds_inds = randperm(noise_count,floor(noise_count*par.noise_level));
        noise_inds = noise_inds_all(noise_inds_inds);
    end
    N = zeros(size(X));
    N(noise_inds) = 1; % TODO: add option for non-binary noise
    X = xor(X,N); % TODO: add option for weighted graphs
    
    if ~strcmp(par.type,'')
        cur_size = par.size(1)+par.size(end)*par.is_bipartite;
        switch par.type
            case {'blank','dummy_blank'}
                par.labels = 1:cur_size;
            otherwise
                par.labels = ones(1,cur_size);
        end
    else
        par.size = size(X,1);
    end

end
