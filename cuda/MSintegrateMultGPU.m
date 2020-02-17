function [u] = MSintegrateMultGPU(gpuFunc,u,throughputVmf_v,rotatedMovmf_v,el,randPhase)
    elDim = size(el);
    elDim(end+1:4) = 1;
    u = feval(gpuFunc.MSintegrateMult,u,size(u),randPhase, el, elDim, ...
        throughputVmf_v.mu1, throughputVmf_v.mu2, throughputVmf_v.mu3, throughputVmf_v.c, throughputVmf_v.dim, ...
        rotatedMovmf_v.mu1, rotatedMovmf_v.mu2, rotatedMovmf_v.mu3, rotatedMovmf_v.c, rotatedMovmf_v.alpha, rotatedMovmf_v.dim);
end

