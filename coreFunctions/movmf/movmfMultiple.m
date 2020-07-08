function [movmfMult] = movmfMultiple(movmf1,movmf2,conjMultiple)
% multiple two movmf components
%
% INPUT
% movmf1,movmf2 : mixture of von Mises-Fisher structres
% conjMultiple: active conjegate while multiplying
%
% OUTPUT
% movmf: mixture of vms after multipication, k equals to k1 * k2, and N is
% according to the matrixMultiple flag

if(nargin < 3)
    conjMultiple = true;
end


if(conjMultiple)
%     movmfMult.alpha = movmf1.alpha .* conj(movmf2.alpha);
    movmfMult.c     = movmf1.c      + conj(movmf2.c)    ;
    movmfMult.mu1   = movmf1.mu1    + conj(movmf2.mu1)  ;
    movmfMult.mu2   = movmf1.mu2    + conj(movmf2.mu2)  ;
    movmfMult.mu3   = movmf1.mu3    + conj(movmf2.mu3)  ;
else
%     movmfMult.alpha = movmf1.alpha .* movmf2.alpha      ;
    movmfMult.c     = gather(movmf1.c)      + gather(movmf2.c)          ;
    movmfMult.mu1   = gather(movmf1.mu1)    + gather(movmf2.mu1)        ;
    movmfMult.mu2   = gather(movmf1.mu2)    + gather(movmf2.mu2)        ;
    movmfMult.mu3   = gather(movmf1.mu3)    + gather(movmf2.mu3)        ;
end

% movmfMult.dim = max(movmf1.dim,movmf2.dim);
% movmfMult.dim = size(movmfMult.alpha);
	
end