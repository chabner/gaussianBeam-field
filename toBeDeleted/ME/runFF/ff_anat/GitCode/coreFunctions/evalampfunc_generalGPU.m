function af_v = evalampfunc_generalGPU(cosang, sct_type,ampfunc,dim)
    % Pre-calculate single scattering rotation amplitude
    Nv=length(cosang);

    switch sct_type
       case 1
           % isotropic scattering
           af_v = gpuArray(ones(1,Nv)/sqrt(2^(dim-1)*pi));

       case 2
           % provided amplitude function
           ang = acos(cosang);
           af_v = gpuArray(evalampfunc_amplitude(ampfunc,ang,dim));

       case {3,4}
           % HG scattering
           af_v = sqrt(evaluateHGGPU(cosang, ampfunc, 1, dim)); 
    end
end