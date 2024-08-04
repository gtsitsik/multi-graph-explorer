function X = parafac2full(A,B,C,mtimesx_exists)
% PARAFAC2FULL Compose 3-mode PARAFAC.
%
%    X = PARAFAC2FULL(A,B,C) returns a tensor X which is the composition of
%    a PARAFAC with factor matrices A,B and C.
%
%    X = PARAFAC2FULL(A,B,C,mtimesx_exists) allows avoiding the expensive 
%    exist operation for 'mtimesx'. This is useful when PARAFAC2FULL is 
%    called multiple times.
if nargin==3
    mtimesx_exists = exist('mtimesx','file')==3;
end
if mtimesx_exists 
    X = mtimesx(A.*permute(C,[3 2 1]),B');
else

    X = zeros([size(A,1),size(B,1),size(C,1)]);
    for i =1:size(C,1)
        X(:,:,i) = (A.*C(i,:))*B';
    end
end

% Third way of calculating X_rec using tensor toolbox. Can be very slow.
% X = double(ktensor({A,B,C}));