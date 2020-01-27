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

movmf.mu = movmf.mu(n,:,:);
movmf.alpha = movmf.alpha(n,:);
movmf.c = movmf.c(n,:);
movmf.N = size(movmf.alpha,1);

end