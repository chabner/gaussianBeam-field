function [movmf] = mogToMovmf(mog)
    movmf.mu1 = mog.b1;
    movmf.mu2 = mog.b2;
    movmf.mu3 = mog.a;
    movmf.c = mog.c - mog.a;
    movmf.alpha = mog.alpha;
    movmf.dim = mog.dim;
    
    if(isfield(mog,'dirDim'))
        movmf.dirDim = mog.dirDim;
    end
end

