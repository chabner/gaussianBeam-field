function [movmf] = movmfPickMixtureComponent(movmf,k_idx)
% pick movmf mixture component along k possible mixtures
%
% INPUT
% movmf: mixture of von Mises-Fisher structre
% k_idx: the mixtures to be picked, can be a vector as long as it is less
% then k
%
% OUTPUT
% movmf: the mixtures which picked

movmf.mu = movmf.mu(:,k_idx,:);
movmf.alpha = movmf.alpha(:,k_idx);
movmf.c = movmf.c(:,k_idx);
movmf.k = size(movmf.alpha,2);

end