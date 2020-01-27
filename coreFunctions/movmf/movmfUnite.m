function [movmf] = movmfUnite(movmf1,movmf2)
% unite two movmf to a single vmf
% the union is on n
%
% INPUT
% movmf1,movmf2: mixture of von Mises-Fisher structre with same k and dim
%
% OUTPUT
% movmf: the united movmf

N1 = movmf1.N;
N2 = movmf2.N;

dim = movmf1.dim;
k = movmf1.k;

movmf.mu = zeros(N1+N2,k,dim);
movmf.alpha = zeros(N1+N2,k);
movmf.c = zeros(N1+N2,k);

movmf.mu(1:N1,:,:) = movmf1.mu;
movmf.mu(N1+1:(N1+N2),:,:) = movmf2.mu;

movmf.alpha(1:N1,:) = movmf1.alpha;
movmf.alpha(N1+1:(N1+N2),:) = movmf2.alpha;

movmf.c(1:N1,:) = movmf1.c;
movmf.c(N1+1:(N1+N2),:) = movmf2.c;

movmf.N = N1 + N2;
movmf.k = k;
movmf.dim = dim;

end