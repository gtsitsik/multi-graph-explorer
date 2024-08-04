function param_sweep_explore_demo()

% Describes an artificial graph with 3 view clusters of 3 views each. The
% first view cluster has 3 node clusters and the other 2 view cluters have
% 2 node clusters.
sz= {...
    {60,40,20},
    {100,20},
    {20,100}
    };
sp = [0.94 0.93 0.92];
G = graph_tree_root;
for i = 1:numel(sz)
    G.Children(i) = graph_tree_node;
    cur_ch = G.Children(i);
    cur_ch.slices_num = 3;
    cur_ch.noise_level = 0.01;
    cur_ch.sparsity_level = sp(i);
    for j = 1:numel(sz{i})
        cur_ch.Children(j).type = 'clique';
        cur_ch.Children(j).size = sz{i}{j};
    end
end
create_graph(G);

% % For use with real-wold data the user needs to provide the adjacency
% % tensor X, and the nodes and view labels. If no labels are available, you
% % can assign all nodes and/or views the same label. For example, if we
% % extract the following information from graph G above 
% X = create_graph(G);
% view_labels = G.labels;
% node_labels = {G.Children.labels};
% % then we can generate a "real-world" version of it as follows
% T = graph_tree_root;
% T.Data = X;
% T.labels = view_labels
% for i=1:numel(node_labels)
%     T.Children(i) = graph_tree_node;
%     cur_ch = T.Children(i);
%     cur_ch.labels = node_labels{i};
% end

% Parameter tree based on which all the parameter combinations will be
% generated
params.graph_tree = G;
params.L_type_ind = [1 2]; % 1:Adjacency matrix,2:Normalized Laplacian
params.embedding_method.ComClus.beta = [0.1,0.9];
params.embedding_method.ComClus.rho = [0,0.16];
params.embedding_method.ComClus.thres_inner = [1e-6];
params.embedding_method.Symmetric_Richcom.structure = "true";
params.embedding_method.Symmetric_Richcom.rho = [0,0.16]; 
params.thres = [1e-6];
params.max_iters = [1000];
params.R = 1:5;
params.M = 3;
params.sample = 1:10; 
params.workers = feature('numcores');  
params.column_normalization_type = ["B"];
params.clustering_method.large_inner_prod.thres = [0.5,0.75,0.9];
params.clustering_method.kmeans.replicates = [3];
params.clustering_method.kmeans.clusters_num = "3 2 2";
params.clustering_method.kmeans.row_normalization_type = ["unit"];  
params.clustering_measure = ["NMI","silhouette_equal"];
params.clustered_entity = "nodes";
[~,~,~,~,~,filespath] = multi_run_calculate(params,true,true);

% Assignment of parameters to the preferred visual dimensions 
fixed_params =[];
xaxis = "R";
tile = ["clustering_measure","L_type_ind"];
ln_clr = ["embedding_method","clustering_method"];
aggregate = "sample";
filename = "cur_plot";
[my_fig,plot_params] = multi_run_GUI(filespath,xaxis,tile,ln_clr,aggregate,fixed_params);
savefig(my_fig,filespath+"/"+filename)
save(filespath+"/"+filename,'plot_params','fixed_params','-v7')

% % Timing measurements for each parameter combination in the generated
% % figure
% samples_num = 5;
% save_is_on = true;
% create_durations_figure(filespath,filename,samples_num,save_is_on)
% end
