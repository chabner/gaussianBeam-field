function [mog] = movmfToMog(movmf)
    mog.b1 = movmf.mu1;
    mog.b2 = movmf.mu2;
    mog.a = movmf.mu3;
    mog.c = movmf.c + movmf.mu3;
    mog.alpha = movmf.alpha;
    mog.dim = movmf.dim;
    
    if(isfield(movmf,'dirDim'))
        mog.dirDim = movmf.dirDim;
    end
end

