function [us,u] = ff_single(gpuFunc,us,u,Wl,e_l0,af_ang_vl,x,constCont)
    setConstantMemory(gpuFunc.ff_single,'fastConstCopy',[x;real(constCont);imag(constCont)]);
%     setConstantMemory(gpuFunc.ff_single,'x',x);
%     setConstantMemory(gpuFunc.ff_single,'constCont',constCont);
    
    [us,u] = feval(gpuFunc.ff_single, us, u, Wl, e_l0, af_ang_vl, ...
        gpuFunc.v.x, gpuFunc.v.y, gpuFunc.v.z, ...
        gpuFunc.dirv.x, gpuFunc.dirv.y, gpuFunc.dirv.z);
end
