function e = evalphaseattAp(x,v,is_ff,sigt,lambda,box_min,box_max, ap_win)
    % complex transmission and attenuation term in multiple scattering
    
    Nv = size(v,2);
    dim = size(v,1);
    if (is_ff)
        pv = x'*v;
        rv = ones(1,Nv);
        [bdv,bx] = cubeDist(x,box_min,box_max,-v);
    else
        d = repmat(x,1,Nv)-v;
        nd = sum(d.^2,1).^0.5;
        v = d./repmat(nd,dim,1);
        pv = nd; 
        rv = 1./(nd+0.01*lambda).^((dim-1)/2);
        [bdv,bx]=cubeDist(x,box_min,box_max,-v);
        bdv = min(bdv,nd);  
    end

     e = exp(1i*pi*2/lambda*pv-sigt/2*bdv).*rv;
     apd=size(ap_win,1);
     m=zeros(1,Nv);
     for jj=1:apd
         m=m.*bx(jj,:)<=ap_win(jj,2);
         m=m.*bx(jj,:)>=ap_win(jj,1);
         
     end
     e=e.*m;
     
end