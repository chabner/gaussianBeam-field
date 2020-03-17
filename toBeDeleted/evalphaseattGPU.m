function [e_v0] = evalphaseattGPU(gpuFunc,x)
    setConstantMemory(gpuFunc.evalphaseattGPU,'x',x);
    
    e_v0 = feval(gpuFunc.evalphaseattGPU,gpuFunc.e_v0, ...
        gpuFunc.v.x, gpuFunc.v.y, gpuFunc.v.z, ...
        gpuFunc.dirv.x, gpuFunc.dirv.y, gpuFunc.dirv.z);
end

