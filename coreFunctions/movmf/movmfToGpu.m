function [movmf_gpu] = movmfToGpu(movmf)
    movmf_gpu = movmf;
    
    movmf_gpu.mu1 = gpuArray(movmf_gpu.mu1);
    movmf_gpu.mu2 = gpuArray(movmf_gpu.mu2);
    movmf_gpu.mu3 = gpuArray(movmf_gpu.mu3);
    movmf_gpu.c = gpuArray(movmf_gpu.c);
end

