function [X,labels,slice_structure_id_orig] = parse_dataset(data_ID)

switch data_ID
    case 1
        fid = fopen('.\air-multi-public-dataset\network.txt');
        labels = [];
        slice_structure_id_orig = zeros(1,size(X,3));
        k = 1;
        
        while true
            cur_line = str2num(fgetl(fid));
            str2double(fgetl(fid));
            for i=1:cur_line
                C = str2num(fgetl(fid));
                for j=1:C(2)
                    X(C(1),C(j),k)=1;
                end
            end
            if k==37
                break ;
            else
                fgetl(fid);
                k=k+1;
            end
        end
        
    case 2
        fid = fopen('.\eu-core\email-Eu-core-department-labels.txt');
        
        
        while ~feof(fid)
            C = str2num(fgetl(fid));
            labels{1}(C(1)+1)=C(2);
        end
        slice_structure_id_orig = 1;
        
        fid = fopen('.\eu-core\email-Eu-core.txt');
        while ~feof(fid)
            C = str2num(fgetl(fid));
            X(C(1)+1,C(2)+1)=1;
        end
        
    case 3
        if isfile("conferences_dblp/conf.mat")
            data = load("conferences_dblp/conf.mat");
            X = data.X;
            labels = data.labels;
            slice_structure_id_orig = data.slice_structure_id_orig;
        else
            X=sptensor([1 1 1]);
            file_names = ["AAAI","ICDE","ICDM","ICML","IJCAI","KDD","NIPS","SIGMOD","VLDB"];
            slice_structure_id_orig=zeros(1,numel(file_names));
%             3 2 1 3 3 1 3 2 2
            labels = [];
            IDs = [];
            n = 1;
            for conf_id = 1:numel(file_names)

                text = fileread("conferences_dblp/"+file_names(conf_id)+".json");
                kdd = jsondecode(text);
            
                for i =1:100%:numel(kdd.result.hits.hit)

                    if mod(i,100)==0
                        disp("i="+i+"/"+numel(kdd.result.hits.hit)+" - conf "+conf_id+"IDs num "+size(IDs,1))
                    end
                    if isfield(kdd.result.hits.hit(i).info,"authors")
                        tmp1 = kdd.result.hits.hit(i).info.authors.author;
                        for j = 1:numel(tmp1)
                            cur_id_1 = [];
                            for k = 1:numel(tmp1(j).x_pid)
                                cur_id_1= [cur_id_1 , num2str(double(tmp1(j).x_pid(k))) ];
                            end
                            cur_id_1 = str2num(cur_id_1);
                            for l = j+1:numel(tmp1)
                                cur_id_2 = [];
                                for m = 1:numel(tmp1(l).x_pid)
                                    cur_id_2= [cur_id_2 , num2str(double(tmp1(l).x_pid(m))) ];
                                end
                                cur_id_2 = str2num(cur_id_2);
                                IDs(n,:) = [cur_id_1 cur_id_2 conf_id ];
                                n=n+1;
                            end
                            
                        end
                    end
                end
                
            end
              disp("-1")
            [uniq_IDs,~,IDs_inds] = unique(reshape(IDs(:,[1 2]),[],1));
            disp("0")
%             IDs(:,1)= int64(sum((IDs(:,1)==uniq_IDs').*[1:numel(uniq_IDs)],2)) ;
%             IDs(:,2)= int64(sum((IDs(:,2)==uniq_IDs').*[1:numel(uniq_IDs)],2));
IDs(:,[1 2])=reshape(IDs_inds,[],2);
             disp("1")
            %             X=zeros(numel(uniq_IDs),numel(uniq_IDs));
            size(IDs)
            for i = 1:size(IDs,1)
                if mod(i,1000)==0
                    disp("i="+i+"/"+size(IDs,1))
                end
                if IDs(i,1)>size(X,1) || IDs(i,2)>size(X,2) || IDs(i,3)>size(X,3)
                    X(IDs(i,1),IDs(i,2),IDs(i,3)) =0;
                end
                if IDs(i,1)>size(X,2) || IDs(i,2)>size(X,1)
                    X(IDs(i,2),IDs(i,1),IDs(i,3)) =0;
                end
                
                X(IDs(i,1),IDs(i,2),IDs(i,3))=X(IDs(i,1),IDs(i,2),IDs(i,3))+1;
                X(IDs(i,2),IDs(i,1),IDs(i,3))=X(IDs(i,1),IDs(i,2),IDs(i,3));
            end
              disp("2")
%             X = padarray(X,[(max(size(X))-size(X,1)) (max(size(X))-size(X,2))],0,'post');%check validity
%             X=X+permute(X,[2,1,3]);
            save("conferences_dblp/conf.mat",'X','slice_structure_id_orig','labels');
        end
        
    otherwise
        disp("wrong dataset")
end


% X = padarray(X,[(max(size(X))-size(X,1)) (max(size(X))-size(X,2))],0,'post');
% if size(X,1)>size(X,2)
%     X = cat(2,X,zeros(size(X,1),size(X,1)-size(X,2)));
% else
%     X = cat(1,X,zeros(size(X,2)-size(X,1),size(X,2)));
% end
for i =1:numel(labels)
    if ~isempty(labels{i})
        ss = [2 3 4 5]';
        ss = [ 1 2 3]';
        inds = find(prod(labels{i}~=ss,1));
        X(inds,:)=[];
        X(:,inds)=[];
        labels{i}( inds )=[];
        labels{i} = sum((unique(labels{i})'==labels{i}).*[1:numel(unique(labels{i}))]',1);
    end
end
for i = 1:numel(labels)
    labels{i}=labels{i}';
end
% tmp = 0;
% for i = 1:size(X,3)
%     XX = X(:,:,i);
%     tmp = tmp + sum(XX,1) + sum(XX,2)' - 2*diag(XX)';
% end
%
% tmp = find(tmp==0);
% X(tmp,:,:)=[];
% X(:,tmp,:)=[];
%
% if ~isempty(labels{i})
%     labels{i}(tmp)=[];
% end