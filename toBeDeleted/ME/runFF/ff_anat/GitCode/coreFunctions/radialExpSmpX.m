function [x,px]=radialExpSmpX(box_min,box_max,l,sigt,lambda,lsign)
    % exponential distribution position of first scatter
    dim=size(l,1);
    Nl=size(l,2);
    Ns0=150;
    Ns=round(Ns0/Nl);
    x=[]; px=[];
    
    V=prod(box_max-box_min);
    
    for j=1:Nl
       tNs=0;
       while(tNs==0)
       w=randn(dim,Ns);
       w=w./repmat(sqrt(sum((w.^2))),dim,1);
       tl=l(:,j);
       [tl,pd]=cubeProj(tl,box_min,box_max,w);
       nni=find(~isnan(tl(1,:)));
       tl=tl(:,nni);
       w=w(:,nni);
       pd=pd(:,nni);
       tNs=length(nni);
       end
       bdw=zeros(1,tNs);
       %pd(:)=0;
       for jj=1:tNs
       bdw(jj) = cubeDist(tl(:,jj),box_min,box_max,w(:,jj));
       end
       pw=exp(-bdw*sigt);
       bdw=bdw+pd;
       rmax=exp(-bdw*sigt);
       rmin=exp(-pd*sigt);
       
       %pw=rmin-rmax;
       
       nii=sampleVec(pw/sum(pw),ceil(tNs/2));
       w=w(:,nii);
       rmin=rmin(:,nii);
       rmax=rmax(:,nii);
       bdw=bdw(:,nii);
       pd=pd(:,nii);
       tl=tl(:,nii);
       tNs=length(nii);
       
       
       d=-log(rand(1,tNs).*(rmax-rmin)+rmin)/(sigt);
       if lsign(j)<0
          d=bdw-d;
       end
       x=[x,l(:,j)+repmat(d,dim,1).*w];
       %px=[px,rmin-rmax];
       %px=[px,1./(d*2*pi).^((dim-1)/2).*exp(-sigt*d0)./(rmax).*sigt];
    end   
    
    
    
    Ns0=size(x,2);
    qL=zeros(Ns0,Nl);
    for j=1:Ns0
       qL(j,:) = evalphaseatt(x(:,j),l,0,sigt,lambda,box_min,box_max,lsign);
    end
    q=(prod(abs(qL),2)).^(1/(Nl/2));
    p=mean(abs(qL).^2,2);
    q=(q./p);%.*px(:);
       
    q=q/sum(q);
    cq=cumsum(q);
    r=rand;
    ii=min(find(r<=cq));
    if (isempty(ii))
        keyboard
    end
    if 0%p(ii)<0.0001
        keyboard
    end
    x=x(:,ii);
    px=p(ii)*V;   
       
 
end