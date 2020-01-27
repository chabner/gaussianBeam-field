function [p] = movmfPdf(movmf,x)
% calculate the pdf of movmf in point x
%
% INPUT
% movmf: mixture of von Mises-Fisher structre
% x    : [M,movmf.dim] desired direction to calculate the pdf
%
% OUTPUT
% p: [movmf.N,M] the pdf result

mu = permute(movmf.mu,[1,4,2,3]);     % mu in size of [N,1,k,dim]
x = permute(x,[3,1,4,2]);             % x in size of [1,M,1,dim]
muTimesX = sum(mu.*x,4);              % muTimesX in size of [N,M,k]
alpha = permute(movmf.alpha,[1,3,2]); % alpha in size of [N,1,k]
c = permute(movmf.c,[1,3,2]);         % c in size of [N,1,k]

p = sum(alpha .* exp(muTimesX + c) , 3);

end

