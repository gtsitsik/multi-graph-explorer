function [U,A,B,iters,time,obj_cur] = symmetric_richcom(X,R,M,thres,max_iters,rho,B,print_type,mtimesx_exists)
% rng(3)
tic
% The optimization steps in this implementation are derived as modified
% versions of the optimization steps of ComClus as presented in the 
% original paper that introduced it. Therefore, input arguments are renamed
% to follow the notation of ComClus.
if exist('X','var')
    A=X;
end
if exist('R','var')
    h=R;
end
if exist('M','var')
    k=M;
end
if exist('B','var')
    S_in=B';
end

% TODO: Improve input validation. Particularly, make sure that M and B are
% handled appropriately when both are given as input. A solution could be
% to crop or augment B so that its dimension matches M.
if nargin ==0
    U_true = rand(10,4);
    S_true = [1 1 0 0;
            0 0 1 1]';
    V_true = [1 1 1 0 0 0;
              0 0 0 1 1 1]';
    A = parafac2full(U_true,U_true,V_true*S_true');
end
if nargin<=1
    if exist('S_true')
        h = size(S_true,1);
    else
        h = 10;
    end
end
if nargin<=2
    if exist('S_true')
        k = size(S_true,2);
    else
        k = 3;
    end
end
if nargin<=3
    thres = 1e-6;
end
if nargin<=4
    max_iters = 2000;
end
if nargin<=5
    rho = 0.1;
end
if nargin<=6
    %     S_in=rand(h,k);
    %     S_max = max(S_in,[],1);
    %     S_in = double(S_in==S_max);
    if exist('S_true')
        S_in = S_true;
    else
        S_in = double([1:k]'==randi(k,1,h))';
    end
    B = S_in';
end
if nargin<=7
    print_type = "basic";
end
if nargin<=8
    mtimesx_exists = exist('mtimesx')==3;
end


U = rand(size(A,1),h);
V = rand(size(A,3),k);

A_3 = reshape(permute(A,[3,1,2]),size(A,3),[])';

iters = 0;

obj_cur = obj_eval(A,U,V,S_in,rho,mtimesx_exists);
obj_prev = inf;
obj_cur_ = obj_cur;
obj_prev_ = obj_prev;

print_progress("begin")
while abs(obj_cur-obj_prev)>=thres*obj_prev && iters<max_iters
    iters = iters +1;
    W = S_in*V';

    if mtimesx_exists
        UW_all = U.*permute(W,[3 1 2]);
        denominator = sum(mtimesx(mtimesx(UW_all,U'),UW_all),3);
        numerator = sum(mtimesx(A,UW_all),3);
    else
        numerator = 0;
        denominator = 0;
        for j =1:size(A,3)
            UW = (U.*W(:,j)');
            denominator = denominator + UW*U'*UW;
            numerator = numerator + A(:,:,j)*UW;
        end
    end


    
    U = U.*((4*numerator./(4*denominator+rho)).^(1/4));
    %         print_progress("U")

    W = A_3;

%    Same as S = khatrirao(U,U)*S_in;
    if mtimesx_exists
        S = reshape(permute(mtimesx(U.*permute(S_in,[3,1,2]),U'),[3,1,2]),size(S_in,2),[])';
    else
        tmp = zeros(size(U,1),size(U,1),size(S_in,2));
       
        for i = 1:size(tmp,3)
            tmp(:,:,i) = (U.*S_in(:,i)')*U';
        end
        S = reshape(permute(tmp,[3,1,2]),size(S_in,2),[])';
    end
    
    V = V.*(((2*W'*S)./(2*V*(S'*S)+rho)).^(1/2));
    %         print_progress("V")
    
    
    obj_prev = obj_cur;
    obj_cur = obj_eval(A,U,V,S_in,rho,mtimesx_exists);
    
    
    print_progress("")
end
print_progress("end")

    function print_progress(var)
        if print_type~="nothing"
            obj_prev_ = obj_cur_;
            X_rec_ = parafac2full(U,U,V*S_in');
            obj_cur_ = obj_eval(A,U,V,S_in,rho,mtimesx_exists);
            obj_change_ = (obj_cur_-obj_prev_)/obj_prev_*100;
            msg = "";
            if max(print_type == ["all","basic"])
                msg = var+" iters: " + iters + " - rec error:   "+num2str(norm(A(:)-X_rec_(:))/norm(A(:))*100) + "   " +  obj_change_;
                if print_type == "all"
                    msg = "Symmetric Richcom R:"+h+" M:"+k+" thres:"+thres+" max_iters:"+max_iters +" rho: "+ rho + " | "+msg;
                end
            end
            if (obj_change_>-inf) || var=="begin" || var=="end"|| var==""
                disp(msg);
            end
        end
    end

% The optimization steps in this implementation are derived as modified
% versions of the optimization steps of ComClus as presented in the 
% original paper that introduced it. Therefore, input arguments are renamed
% to follow the notation of ComClus.
A = V;
time = toc;
end

function L_all = obj_eval(A,U,V,S,rho,mtimesx_exists)
X_rec = parafac2full(U,U,V*S',mtimesx_exists);
L_rec = norm(A(:)-X_rec(:))^2;
% L_R = norm(W-S*V','fro')^2;
L_reg = norm(V(:),1)+norm(U(:),1);

L_all = L_rec  + rho*L_reg;
end
