function [u] = MSintegrateMultGPU(gpuFunc,u,throughputVmf_v,rotatedMovmf_v,el,randPhase)
    setConstantMemory(gpuFunc.MSintegrateMult,'randPhase',2*pi*randPhase);

    u = feval(gpuFunc.MSintegrateMult,u, gpuArray(el), ...
        throughputVmf_v.mu1, throughputVmf_v.mu2, throughputVmf_v.mu3, throughputVmf_v.c, ...
        rotatedMovmf_v.mu1, rotatedMovmf_v.mu2, rotatedMovmf_v.mu3, rotatedMovmf_v.c);
end

