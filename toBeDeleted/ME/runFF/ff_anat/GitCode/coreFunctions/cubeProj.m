function [xd,d]=cubeProj(x,vmin,vmax,v)

    thr=0.1^7;
    bd = inf(size(v));
    for j=1:size(vmin,1)
        bd(j,:)=(x(j)-vmax(j))./max(abs(v(j,:)),thr).*(v(j,:)<=0)+(-x(j)+vmin(j))./max(abs(v(j,:)),thr).*(v(j,:)>0);
    end
  
    d = max(bd);

    xd=x*ones(1,size(v,2));
    for j=1:size(v,2)
    if d(j)>=0
        ep=0.1^4*max(abs(vmax-vmin));

        xd(:,j)=x+(d(j)+ep)*v(:,j);
        
         if ~prod((xd(:,j)<=vmax).*(xd(:,j)>=vmin))
             xd(:,j)=nan;
         end
    else
        if prod((x<=vmax).*(x>=vmin))
        xd(:,j)=x;
        else
            xd(:,j)=nan;
        end
    end
    end

    d=max(d,0);
    
end