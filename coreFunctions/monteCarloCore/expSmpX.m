function [x,px]=expSmpX(box_min,box_max,lmean,sigs,smpNum)
    if(nargin < 5)
        smpNum = 1;
    end

    % exponential distribution position of first scatter
    dim=size(lmean,1);
    x=zeros(dim,1,smpNum);   
    lz=abs(lmean(end));
    sigs=sigs/lz;
    dmax=(box_max(end)-box_min(end));
    rmax=1-exp(-dmax*sigs);
    d=-log(-rand(1,1,smpNum).*rmax+1)/(sigs);
    x(1:end-1,1,:) = rand([dim-1,smpNum]).*(box_max(1:end-1)-box_min(1:end-1))+box_min(1:end-1);
    x(end,1,:)=(d+box_min(end)).*(lmean(end)>0)+(-d+box_max(end)).*(lmean(end)<=0);

    px = exp(-sigs*d)/(rmax)*sigs*dmax;
end