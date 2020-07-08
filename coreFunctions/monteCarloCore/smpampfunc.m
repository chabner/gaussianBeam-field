function ang=smpampfunc(ampfunc,dim,smpNum)
    if(nargin < 3)
        smpNum = 1;
    end
    
    % sample direction by provided amplitude function
    N=length(ampfunc.sampleIcdf);
    stp=1/(N-1);
    a=rand(1,1,smpNum);
    i=round(a/stp)+1;
    ang=ampfunc.sampleIcdf(i) * pi/N;
    ang = permute(ang , [1,6,3,4,5,2]);
    if(dim == 2)
        ang = ang * 2;
    end
end
