function [U,A,B,W,iters,time,obj_cur] = comclus(X,R,M,thres,max_iters,beta,rho,thres_inner,print_type,mtimesx_exists)


tic;
% Rename variables for consistency
if exist('X')==1
    A=X;
end
if exist('R')==1
    h=R;
end
if exist('M')==1
    k=M;
end



if nargin ==0
    U_true = rand(15,4);
    S_true = [1 1 0 0;
        0 0 1 1]';
    V_true = [1 1 1 0 0 0;
        0 0 0 1 1 1]';
    A = parafac2full(U_true,U_true,V_true*S_true');
end
if nargin<=1
    h = size(S_true,1);
end
if nargin<=2
    k=size(S_true,2);
end
if nargin<=3
    thres = 1e-6;
end
if nargin<=4
    max_iters=2000;
end
if nargin<=5
    beta = 1;
end
if nargin<=6
    rho = 0.1;
end
if nargin<=7
    thres_inner = thres;
end
if nargin<=8
    print_type ="nothing";
end
if nargin<=9
    mtimesx_exists = exist('mtimesx')==3;
end

U = rand(size(A,1),h);
V = rand(size(A,3),k);
S = rand(h,k);
W = rand(h,size(A,3));

iters = 0;
iters_inner =0;

obj_cur = obj_eval(A,U,V,S,W,beta,rho,mtimesx_exists);
obj_prev = inf;
obj_cur_ = obj_cur;
obj_prev_ = obj_prev;


print_progress("begin")
while abs(obj_cur-obj_prev)>=thres*obj_prev && iters<max_iters
    % it = 0
    iters = iters + 1;
    iters_inner = 0;
    
    while abs(obj_cur-obj_prev)>=thres_inner*obj_prev && iters_inner<max_iters
        
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
        
        % Third way of calculating numerator and denominator.
        
        %         UW_all = U.*permute(W,[3 1 2]);
        %         UW_mat_1=reshape(UW_all,size(UW_all,1),[]);
        %         UW_mat_2=reshape(permute(UW_all,[2 1 3]),size(UW_all,2),[]);
        %         denominator= UW_mat_1*kron(eye(size(A,3)),U')*UW_mat_2';
        %         numerator = reshape(A,size(A,1),[])*UW_mat_2';
        
        
        U = U.*((4*numerator./(4*denominator+rho)).^(1/4));
        %         print_progress("U")
        
        V = V.*(((2*beta*W'*S)./(2*beta*V*(S'*S)+rho)).^(1/2));
        %         print_progress("V")
        
        S = S.*(((2*beta*W*V)./(2*beta*S*(V'*V)+rho)).^(1/2));
        %         print_progress("S")
        
        obj_prev = obj_cur;
        obj_cur = obj_eval(A,U,V,S,W,beta,rho,mtimesx_exists);
        iters_inner = iters_inner +1;
    end
    
    print_progress("U,V,S")
    
    
    
    UtU =U'*U;
    if mtimesx_exists
        tmp2 = mtimesx(mtimesx(U',A),U);
        for p = 1:h
            tmp = mtimesx((UtU.*permute(W,[3,1,2])),UtU);
            z2 = UtU(p,p)^2 + beta;
            z1 = permute(tmp(p,p,:),[1,3,2]) - permute(tmp2(p,p,:),[1,3,2]) + beta*W(p,:) - beta*S(p,:)*V';
            W(p,:) = max([W(p,:)-z1/z2 ;zeros(1,size(W,2))],[],1);
        end
    else
        for i = 1:size(A,3)
            tmp2 = U'*A(:,:,i)*U;
            for p = 1:h
                tmp = (UtU.*W(:,i)')*UtU;
                z2 = UtU(p,p)^2 + beta;
                z1 = tmp(p,p) - tmp2(p,p) + beta*W(p,i) - beta*S(p,:)*V(i,:)';
                W(p,i) = max([W(p,i)-z1/z2 , 0]);
            end
        end
    end
    
    print_progress("W    ")
    
    obj_prev = obj_cur;
    obj_cur = obj_eval(A,U,V,S,W,beta,rho,mtimesx_exists);
end
print_progress("end")

    function print_progress(var)
        if print_type~="nothing"
            obj_prev_ = obj_cur_;
            X_rec_ = parafac2full(U,U,W',mtimesx_exists);
            obj_cur_ = obj_eval(A,U,V,S,W,beta,rho,mtimesx_exists);
            obj_change_ = (obj_cur_-obj_prev_)/obj_prev_*100;
            msg="";
            if max(print_type == ["all","basic"])
                msg = var+" iters: " + iters + " iters_inner: " + iters_inner  + " - rec error:   "+num2str(norm(A(:)-X_rec_(:))/norm(A(:))*100) + "   " +  obj_change_;
                if print_type == "all"
                    msg = "ComClus R:"+h+" M:"+k+" thres:"+thres+" max_iters:"+max_iters+" beta: " + beta +" rho: "+ rho + " thres_inner:"+ thres_inner + " | "+msg;
                end
            end
            if (obj_change_>-inf) || var=="begin" || var=="end"|| var==""
                disp(msg);
            end
        end
    end

% Rename variables for consistency
A=V;
B=S';
time = toc;
end

function L_all = obj_eval(A,U,V,S,W,beta,rho,mtimesx_exists)
X_rec = parafac2full(U,U,W',mtimesx_exists);
L_rec = norm(A(:)-X_rec(:))^2;
L_R = norm(W-S*V','fro')^2;
L_reg = norm(V(:),1)+norm(U(:),1)+norm(S(:),1);

L_all = L_rec + beta*L_R + rho*L_reg;
end