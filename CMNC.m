function [U,A,B,iters,time,obj_cur] = CMNC(X,R,M,thres,max_iters,delta,B,print_type,mtimesx_exists)
tic
[g,p_sd,p_gn,tmp22,tmp4,J_U_aux,JJ,tmp5] = deal([]);
if ~exist('X','var') || isempty(X)
    U_true = rand(30,4);
    S_true = [1 1 0 0
        0 0 1 1]';
    V_true = [1 1 1 0 0 0
        0 0 0 1 1 1]';
    X = parafac2full(U_true,U_true,V_true*S_true');
end
if ~exist('R','var') || isempty(R)
    if exist('S_true')
        R = size(S_true,1);
    else
        R = 10;
    end
end
if ~exist('M','var') || isempty(M)
    if exist('S_true')
        M = size(S_true,2);
    else
        M = 3;
    end
end
if ~exist('thres','var') || isempty(thres)
    thres = 1e-6;
end
if ~exist('max_iters','var') || isempty(max_iters)
    max_iters = 2000;
end

if ~exist('delta','var') || isempty(delta)
    delta = 1;
end
if ~exist('B','var') || isempty(B)
    %     S_in=rand(R,M);
    %     S_max = max(S_in,[],1);
    %     S_in = double(S_in==S_max);
    if exist('S_true')
        S_in = S_true;
    else
        S_in = double([1:M]'==randi(M,1,R))';
    end
    B = S_in';
end
if ~exist('print_type','var') || isempty(print_type)
    print_type = "basic";
end
if ~exist('mtimesx_exists','var') || isempty(mtimesx_exists)
    mtimesx_exists = exist('mtimesx')==3;
end



iters = 0;
U_aux = randn(size(X,1),R);
A_aux = randn(size(X,3),M);
JJ = zeros(numel(X),numel(U_aux)+numel(A_aux));
numel_U_aux = 1:numel(U_aux);
U = U_aux.^2;
A = normalize_fibers(A_aux.^2,2);
[obj_cur,r_all] = obj_eval(X,U,A,B,mtimesx_exists);
obj_prev = inf;
obj_cur_ = obj_cur;
obj_prev_ = obj_prev;
steps_need_update = true;
print_progress("begin")
while (abs(obj_cur-obj_prev)>=thres*obj_prev || (~steps_need_update && abs(delta-delta_prev)>=thres*delta_prev)) && iters<max_iters
    iters = iters+1;
    if steps_need_update 
        J = jacobian_eval(X,U_aux,A_aux,U,A,B,r_all);
    end
    p = dogleg(J,r_all,delta,steps_need_update);
    U_aux_new = U_aux + reshape(p(1:numel(U_aux)),size(U_aux'))';
    A_aux_new = A_aux + reshape(p(numel(U_aux)+1:end),size(A_aux'))';
    U_new = U_aux_new.^2;
    A_new = normalize_fibers(A_aux_new.^2,2);
    m_cur = 0.5*norm(J*p+r_all)^2;
    obj_prev = obj_cur;
    r_all_prev = r_all;
    [obj_cur,r_all] = obj_eval(X,U_new,A_new,B,mtimesx_exists);
    gamma = (obj_prev-obj_cur)/(obj_prev-m_cur);

    delta_prev = delta;
    if gamma>0.75 && abs(norm(p)-delta)<32*eps
        delta = 2*delta;
    elseif gamma<0.25
        delta = delta/4;
    end
%     %%%%%%%%%%%%
%         disp("norm(p):                           "+norm(p))
%         disp("delta_prev:                        "+delta_prev)
%         disp("delta:                             "+delta)
%         disp("gamma:                             "+gamma)
%         disp("(obj_prev-obj_cur)/abs(obj_prev)%: "+(obj_prev-obj_cur)/abs(obj_prev)*100)
%         disp("(obj_prev-m_cur)/abs(obj_prev)%:   "+(obj_prev-m_cur)/abs(obj_prev)*100)
%     %%%%%%%%%%
    if gamma>0
        U_aux = U_aux_new;
        A_aux = A_aux_new;
        U = U_new;
        A = A_new;
        steps_need_update = true;
    else
        obj_cur = obj_prev;
        r_all = r_all_prev;
        steps_need_update = false;
    end
    print_progress("U,A")
end
print_progress("end")
time = toc;

    function print_progress(var)
        if print_type~="nothing"
            obj_prev_ = obj_cur_;
            X_rec_ = parafac2full(U,U,A*B,mtimesx_exists);
            obj_cur_ =  obj_eval(X,U,A,B,mtimesx_exists);
            obj_change_ = (obj_cur_-obj_prev_)/obj_prev_*100;

            msg="";
            if max(print_type == ["all","basic"])
                msg = var+" iters: " + iters + " - rec error:   "+num2str(norm(X(:)-X_rec_(:))/norm(X(:))*100) + "   " +  obj_change_;
                if print_type == "all"
                    msg = "CMNC"+" R:"+R+" M:"+M+" thres:"+thres+" max_iters:"+max_iters +" | "+msg;
                end
            end
            if (obj_change_>-inf) || var=="begin" || var=="end"|| var==""
                disp(msg);
                if var=="end"
                    disp("---------------------------------------------------")
                end
            end
        end
    end

    function p_opt = dogleg(J,r_all,delta,steps_need_update) 

        if steps_need_update
            g = J'*r_all; % Gradient
            p_sd = -(g'*g)/norm(J*g)^2*g; % Optimal steepest descend step
            p_gn = -lsqminnorm(sparse(J),r_all); % Gauss-Newton step
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % p_gn = -pinv(J)*r_all;
        % norm(p_gn-p_gn_real)/norm(p_gn_real)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if norm(p_gn)<=delta
            p_opt = p_gn;
        elseif norm(p_sd)>=delta
            p_opt = p_sd/norm(p_sd)*delta;
        else % solving ||p_sd+beta*(p_gn-p_sd)||=delta  
            p_gnsd = p_gn-p_sd;
            a = p_gnsd'*p_gnsd;
            b = 2*p_gnsd'*p_sd;
            c = p_sd'*p_sd-delta^2;
            D = sqrt(b^2-4*a*c);
            beta = (-b+D)/(2*a); %TODO: check validity of this solution
            p_opt = p_sd+beta*(p_gnsd);
        end
        %%%%%%%%%%%%%
%         pl1=[];
%         pl2=[];
%         pl3=[];
%         beta_all = 0:0.01:2;
%         for beta = beta_all
%             if beta<=1
%                 p_test = beta*p_sd;
%             else
%                 p_test = p_sd+(beta-1)*(p_gn-p_sd);
%             end
%             pl1(end+1)= norm(p_test);
%             pl2(end+1) = 0.5*norm(J*p_test+r_all)^2;
%             U_aux_new = U_aux + reshape(p_test(1:numel(U_aux)),size(U_aux'))';
%             A_aux_new = A_aux + reshape(p_test(numel(U_aux)+1:end),size(A_aux'))';
%             U_new = U_aux_new.^2;
%             A_new = normalize_fibers(A_aux_new.^2,2);
%             pl3(end+1) = obj_eval(X,U_new,A_new,B,mtimesx_exists);
%         end
%         plot(beta_all(2:end),sign(diff(pl1/norm(pl1))))
%         hold on
%         plot(beta_all(2:end),0.9*sign(diff(pl2/norm(pl2))))
%         plot(beta_all(2:end),0.8*sign(diff(pl3/norm(pl3))))
%         ylim([-1.2 1.2])
%         hold off
%         pause
        %%%%%%%%%%%%%
    end

    function J = jacobian_eval(X,U_aux,A_aux,U,A,B,r_all)
        tmp = (A*B)';
        tmp = U.*reshape(tmp,[1 size(tmp)]);
        for i = 0:size(U,1)-1
            tmp3 = tmp;
            tmp3(i+1,:,:) = tmp3(i+1,:,:)*2; 
            tmp22(:,i+1,:,i*size(U,2)+(1:size(U,2))) = permute(tmp3,[1,3,2]);
        end

        % Jacobian of the error vector r_all with respect to U
        J_U = reshape(tmp22,numel(r_all),numel(U)); 
        
        % Jacobian of the error vector r_all with respect to U_aux
        J_U_aux = J_U.*reshape(2*U_aux',1,[]);
        JJ(:,1:numel(U_aux)) = J_U_aux; 
        tmp = reshape(U,[size(U,1) 1 size(U,2)]).*reshape(U,[1 size(U)]);
        tmp = reshape(reshape(tmp,[],size(tmp,3))*B',[size(tmp,1,2) size(B,1)]);
        for i = 0:size(A,1)-1
            tmp4(:,:,i+1,i*size(A,2)+(1:size(A,2))) = tmp;
        end

        % Jacobian of the error vector r_all with respect to A
        J_A = reshape(tmp4,numel(r_all),numel(A));

        % Jacobian of the error vector r_all with respect to A_aux
        A_aux_2 = A_aux.*A_aux;
        for i = 0:size(A,1)-1
            tmp = A_aux_2(i+1,:);
            tmp2 = -tmp'*(tmp/norm(tmp)^3);
            tmp2(1:size(tmp2,1)+1:end) = 1/norm(tmp)+diag(tmp2);
            tmp = i*size(A,2)+(1:size(A,2)); 
            tmp5 = J_A(:,tmp)*tmp2.*(2*A_aux(i+1,:));
            JJ(:,numel(U_aux)+tmp) = tmp5; 
        end

        J = JJ; 
    end

end

function [L_all,r_all] = obj_eval(X,U,A,B,mtimesx_exists)
X_rec = parafac2full(U,U,A*B,mtimesx_exists);
r_all = X_rec(:)-X(:);
L_all = 0.5*norm(r_all)^2;
end
