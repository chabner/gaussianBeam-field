function [u] = ff_ms(gpuFunc,u,e_l0_ms,x,w,constCont)
    setConstantMemory(gpuFunc.ff_ms,'fastConstCopy',[x;w;real(constCont);imag(constCont)]);
%     setConstantMemory(gpuFunc.ff_ms,'x',x);
%     setConstantMemory(gpuFunc.ff_ms,'w',w);
%     setConstantMemory(gpuFunc.ff_ms,'constCont',constCont);
    
    u = feval(gpuFunc.ff_ms,u,e_l0_ms, ...
        gpuFunc.v.x, gpuFunc.v.y, gpuFunc.v.z, ...
        gpuFunc.dirv.x, gpuFunc.dirv.y, gpuFunc.dirv.z);
end

