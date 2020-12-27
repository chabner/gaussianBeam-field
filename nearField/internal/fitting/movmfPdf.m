function [p] = movmfPdf(movmf,x)
% calculate the pdf of movmf in point x
%
% INPUT
% movmf: mixture of von Mises-Fisher structre
% x    : [M,3] desired direction to calculate the pdf
%
% OUTPUT
% p: [movmf.N,M] the pdf result

Ndim = numel(movmf.dim);
dim = size(x,2);

if(dim == 3)
    x_1 = x(:,1); x_1 = x_1(:); x_1 = permute(x_1,[2:1:(Ndim+1),1]);
    x_2 = x(:,2); x_2 = x_2(:); x_2 = permute(x_2,[2:1:(Ndim+1),1]);
    x_3 = x(:,3); x_3 = x_3(:); x_3 = permute(x_3,[2:1:(Ndim+1),1]);

    muTimesX = movmf.mu1.*x_1 + movmf.mu2.*x_2 + movmf.mu3.*x_3;
end

if(dim == 2)
    x_1 = x(:,1); x_1 = x_1(:); x_1 = permute(x_1,[2:1:(Ndim+1),1]);
    x_2 = x(:,2); x_2 = x_2(:); x_2 = permute(x_2,[2:1:(Ndim+1),1]);

    muTimesX = movmf.mu1.*x_1 + movmf.mu2.*x_2;
end

p = sum(movmf.alpha .* exp(muTimesX + movmf.c) , 1);

end

