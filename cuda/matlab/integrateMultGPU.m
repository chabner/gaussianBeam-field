function [us] = integrateMultGPU(gpuFunc,us,conv_vmf_l0,throughputVmf_v,randPhase)
    setConstantMemory(gpuFunc.integrateMult,'randPhase',2*pi*randPhase);
    
    us = feval(gpuFunc.integrateMult,us, ...
        conv_vmf_l0.mu1, conv_vmf_l0.mu2, conv_vmf_l0.mu3, conv_vmf_l0.c, ...
        throughputVmf_v.mu1, throughputVmf_v.mu2, throughputVmf_v.mu3, throughputVmf_v.c);
end

