function tests = create_graph_test
tests = functiontests(localfunctions);
end


function test_correct_sparsity_symmetric(testCase)
par = testCase.TestData.par;
for i = 1:numel(par.Children)
    par.Children(i).is_symmetric = true;
    par.Children(i).noise_level = 0.0;
end
X = create_graph(par);
for i = 1:numel(par.Children)
    sparsity_levels(i) = par.Children(i).sparsity_level;
    par.Children(i).sparsity_level = 0.0;
end
X_2 = create_graph(par);
for i = 1:numel(par.Children)
    verifyEqual(testCase,round(1-nnz(X(:,:,par.labels==i))/nnz(X_2(:,:,par.labels==i)),2),sparsity_levels(i))
end
end


function test_correct_sparsity_asymmetric(testCase)
par = testCase.TestData.par;
for i = 1:numel(par.Children)
    par.Children(i).is_symmetric = false;
    par.Children(i).noise_level = 0.0;
end
X = create_graph(par);
for i = 1:numel(par.Children)
    sparsity_levels(i) = par.Children(i).sparsity_level;
    par.Children(i).sparsity_level = 0.0;
end
X_2 = create_graph(par);
for i = 1:numel(par.Children)
    verifyEqual(testCase,round(1-nnz(X(:,:,par.labels==i))/nnz(X_2(:,:,par.labels==i)),2),sparsity_levels(i))
end
end


function test_correct_noise_symmetric(testCase)
par = testCase.TestData.par;
for i = 1:numel(par.Children)
    par.Children(i).is_symmetric = true;
    par.Children(i).sparsity_level = 0.0;
end
X = create_graph(par);
for i = 1:numel(par.Children)
    noise_levels(i) = par.Children(i).noise_level;
    par.Children(i).noise_level = 0.0;
end
X_2 = create_graph(par);
for i = 1:numel(par.Children)
    dim_3_inds = par.labels==i;
    verifyEqual(testCase,round(nnz(X_2(:,:,dim_3_inds)~=X(:,:,dim_3_inds))/(numel(X_2(:,:,dim_3_inds))-numel(dim_3_inds)*size(X_2,1)),2),noise_levels(i))
end
end



function test_correct_noise_asymmetric(testCase)
par = testCase.TestData.par;
for i = 1:numel(par.Children)
    par.Children(i).is_symmetric = false;
    par.Children(i).sparsity_level = 0.0;
end
X = create_graph(par);
for i = 1:numel(par.Children)
    noise_levels(i) = par.Children(i).noise_level;
    par.Children(i).noise_level = 0.0;
end
X_2 = create_graph(par);
for i = 1:numel(par.Children)
    dim_3_inds = par.labels==i;
    cur_X = X(:,:,par.labels==i);
    cur_X_2 = X_2(:,:,par.labels==i);
    verifyEqual(testCase,round(nnz(cur_X_2~=cur_X)/(numel(cur_X_2)-numel(dim_3_inds)*size(cur_X_2,1)),2),noise_levels(i))
end
end

function test_empty_diagonal_symmetric(testCase)
par = testCase.TestData.par;
X = create_graph(par);
for i = 1:size(X,3)
    verifyEqual(testCase,diag(X(:,:,i)),zeros(size(X(:,:,i),1),1))
end
end

function test_empty_diagonal_asymmetric(testCase)
par = testCase.TestData.par;
for i = 1:numel(par.Children)
    par.Children(i).is_symmetric = false;
end

X = create_graph(par);
for i = 1:size(X,3)
    verifyEqual(testCase,diag(X(:,:,i)),zeros(size(X(:,:,i),1),1))
end
end

function test_repeated_use_of_create_graph(testCase)
par = testCase.TestData.par;
par2 = par.deep_copy;
create_graph(par);
create_graph(par2);
verifyTrue(testCase,par.is_equivalent_to(par2))
create_graph(par);
create_graph(par);
create_graph(par);
create_graph(par);
verifyTrue(testCase,par.is_equivalent_to(par2))

end

function setup(testCase)
sizes= [30 20 10;50 10 0]*5;
par = graph_tree_root;
par.Children(1).is_symmetric = true;
par.Children(1).slices_num = 3;
par.Children(1).noise_level = 0.29;
par.Children(1).sparsity_level = 0.55;
par.Children(1).Children(1).size = sizes(1,1);
par.Children(1).Children(1).type = 'clique';
par.Children(1).Children(2).size = sizes(1,2);
par.Children(1).Children(2).type = 'clique';
par.Children(1).Children(3).size = sizes(1,3);
par.Children(1).Children(3).type = 'clique';

par.Children(2).is_symmetric = true;
par.Children(2).slices_num = 1;
par.Children(2).noise_level = 0.33;
par.Children(2).sparsity_level = 0.35;
par.Children(2).Children(1).size = sizes(2,1);
par.Children(2).Children(1).type = 'clique';
par.Children(2).Children(2).size = sizes(2,2);
par.Children(2).Children(2).type = 'clique';
par.Children(2).Children(3).size = sizes(2,3);
par.Children(2).Children(3).type = 'clique';

testCase.TestData.par = par;
end
