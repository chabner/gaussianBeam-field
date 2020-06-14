function [movmf_united] = movmfUnite(movmf_1,movmf_2)
    movmf_united = movmf_1;
    
    movmf_united.mu1 = cat(5,movmf_1.mu1,movmf_2.mu1);
    movmf_united.mu2 = cat(5,movmf_1.mu2,movmf_2.mu2);
    movmf_united.mu3 = cat(5,movmf_1.mu3,movmf_2.mu3);
    movmf_united.c = cat(5,movmf_1.c,movmf_2.c);
    movmf_united.alpha = cat(5,movmf_1.alpha,movmf_2.alpha);
    movmf_united.dir = cat(5,movmf_1.dir,movmf_2.dir);
    
    movmf_united.dim(5) = 2;
end
