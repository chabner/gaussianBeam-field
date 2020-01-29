function af_v = evalampfunc_general(cosang, sct_type,ampfunc,dim)
    % Pre-calculate single scattering rotation amplitude
    Nv=length(cosang);

    switch sct_type
       case 1
           % isotropic scattering
           af_v = ones(1,Nv)/sqrt(2^(dim-1)*pi);

       case 2
           % provided amplitude function
           ang = acos(cosang);
           af_v = evalampfunc_amplitude(ampfunc,ang,dim);

       case {3,4}
           % HG scattering
           af_v = sqrt(evaluateHG(cosang, ampfunc.g, 1, dim)); 
    end
end