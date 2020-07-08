function [movmf_united] = movmfUnite(movmf_1,movmf_2,dimCat)

    if(nargin <3)
        dimCat = 5;
    end
    
    movmf_united = movmf_1;
    
    movmf_united.mu1 = cat(dimCat,movmf_1.mu1,movmf_2.mu1);
    movmf_united.mu2 = cat(dimCat,movmf_1.mu2,movmf_2.mu2);
    movmf_united.mu3 = cat(dimCat,movmf_1.mu3,movmf_2.mu3);
    movmf_united.c = cat(dimCat,movmf_1.c,movmf_2.c);
    
    movmf_united.dir1 = cat(dimCat,movmf_1.dir1,movmf_2.dir1);
    movmf_united.dir2 = cat(dimCat,movmf_1.dir2,movmf_2.dir2);
    movmf_united.dir3 = cat(dimCat,movmf_1.dir3,movmf_2.dir3);
    
%     movmf_united.alpha = cat(dimCat,movmf_1.alpha,movmf_2.alpha);
%     movmf_united.dir = cat(dimCat,movmf_1.dir,movmf_2.dir);
%     
%     movmf_united.dim(dimCat) = 2;
end
