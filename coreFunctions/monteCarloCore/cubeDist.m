function [d,xd] = cubeDist(x,vmin,vmax,v)
    % distance from v to x inside the box

%     thr = 0.1^7;
    v_deminator = abs(v);
%     deminator(deminator < thr) = thr;
    x_vmin = x-vmin;
    x_vmax = -x+vmax;
    bd = (x_vmin./v_deminator).*(v<=0) + (x_vmax./v_deminator).*(v>0);

    d=min(bd);
    
    if nargout>1
        ep = 0.1^4*max(abs(vmax-vmin));
        dim=size(v,1);
        
        xd = x+repmat((d+ep),dim,1).*v;
    end
    
end