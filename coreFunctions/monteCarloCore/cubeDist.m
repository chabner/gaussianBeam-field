function [d,xd] = cubeDist(x,vmin,vmax,v)
    % distance from v to x inside the box

    thr = 0.1^7 * ones(size(vmin,1),1,class(v));
    deminator = max(abs([v,thr]),[],2);
    bd = (x-vmin)./deminator.*(v<=0) + (-x+vmax)./deminator.*(v>0);

    d=min(bd);
    
    if nargout>1
        ep = 0.1^4*max(abs(vmax-vmin));
        dim=size(v,1);
        
        xd = x+repmat((d+ep),dim,1).*v;
    end
    
end