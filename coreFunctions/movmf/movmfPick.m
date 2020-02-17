function [movmf] = movmfPick(movmf,n)
% pick movmf mixtures among the N possible mixtrues
%
% INPUT
% movmf: mixture of von Mises-Fisher structre
% n    : the mixtures to be picked, can be a vector as long as it is less
% then N
%
% OUTPUT
% movmf: the mixtures which picked

movmf.mu1 = movmf.mu1(:,n);
movmf.mu2 = movmf.mu2(:,n);
movmf.mu3 = movmf.mu3(:,n);
movmf.alpha = movmf.alpha(:,n);
movmf.c = movmf.c(:,n);
movmf.dim = size(movmf.alpha);
movmf.dim(end+1:4) = 1;

end