function fig = create_durations_figure(filespath,filename,samples_num,save_is_on)
if ~exist('filespath','var') || ...
        ( ~isstring(filespath) && ~ischar(filespath)) || ...
        ~exist(filespath,'dir')
    error('filespath is not valid')
end
if ~exist('filename','var') || ...
        (~isstring(filename) && ~ischar(filename)) 
    error('filename is not valid')
elseif exist(filespath+"/"+filename+".fig",'file')~=2 
    error(filename+".fig does not exist")
elseif exist(filespath+"/"+filename+".mat",'file')~=2
    error(filename+".mat does not exist")
end

if ~exist('samples_num','var') || isempty(samples_num)
    samples_num=10;
end
if ~exist('save_is_on','var') || isempty(save_is_on)
    save_is_on = false;
end
fig = openfig(filespath+"/"+filename+".fig");
linkaxes(findobj(fig.Children.Children,'Type','Axes'),'off')
load(filespath+"/"+filename+".mat")
% delete(findall(gcf,'Type','Patch'))
my_tiledlayout = fig.Children;
for i = 1:size(plot_params,1)
    ax = nexttile(my_tiledlayout,i);
    ln_all = flip(ax.Children);
    ptch_all = flip(findall(ax,'Type','Patch'));
    ylabel("Duration (sec)")
    for j = 1:size(plot_params,2)
        for k = 1:size(plot_params,3)
            disp("=== Tile:"+i+" | Line: "+j+" | Point: "+k+" ===")
            cur_params = plot_params{i,j,k}.params;
            cur_params.sample = 1:samples_num;
            [~,~,duration,duration_sum] = multi_run_calculate(cur_params,false,false,"time");
%             ln_all(j).YData(k) = duration/samples_num;
            prc = prctile(duration_sum,[25 50 75]);
            ln_all(j).YData(k) = prc(2);
            ptch_all(j).YData(k) = prc(3);
            ptch_all(j).YData(end-k+1) = prc(1);
            axis tight
            drawnow
        end
    end
end
linkaxes(findobj(fig.Children.Children,'Type','Axes'))
if save_is_on
    savefig(fig,filespath+"/"+filename+"_durations_"+samples_num+"samples.fig")
end
end
