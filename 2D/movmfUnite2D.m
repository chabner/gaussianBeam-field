function [movmf_united] = movmfUnite2D(movmf_1,movmf_2)
    movmf_united = movmf_1;
    
    movmf_united.mu1 = cat(movmf_united.focalDim,movmf_1.mu1,movmf_2.mu1);
    movmf_united.mu2 = cat(movmf_united.focalDim,movmf_1.mu2,movmf_2.mu2);
    movmf_united.c = cat(movmf_united.focalDim,movmf_1.c,movmf_2.c);
    movmf_united.alpha = cat(movmf_united.focalDim,movmf_1.alpha,movmf_2.alpha);
    
    movmf_united.dim(movmf_united.focalDim) = ...
        movmf_1.dim(movmf_united.focalDim) + movmf_2.dim(movmf_united.focalDim);
end

