function w = smpampfunc_general(ow, sct_type,ampfunc,smpNum)
    % sample new scatteing direction
    
    if(nargin == 3)
        smpNum = 1;
    end

    dim = size(ow,1);
    switch sct_type
        case 1
        w = randn(dim,1);
        w = w/norm(w);
        
        case 2
        a = smpampfunc(ampfunc,dim,smpNum);
        cosa = cos(a);

        case 3
        cosa = sampleHG(ampfunc,dim,smpNum);
        
        case 4
        cosa = sampleHG(ampfunc,dim,smpNum);
        
        
    end

    if sct_type>1
        w = rotateBytheta(ow,cosa,smpNum);
    end

end