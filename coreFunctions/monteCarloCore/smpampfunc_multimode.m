function [w,pw] = smpampfunc_multimode(ow, sct_type,ampfunc)
    % sample new scatteing direction

    
    Ns=200;
    %Ns=10000;
    Ns=10;
    dim = size(ow,1);
    switch sct_type
         
        case 2
            a=zers(1,Ns);
            for j=1:Ns
              a(j) = smpampfunc(ampfunc);
            end
            cosa = cos(a);

        case {3,4}
            cosa=zeros(1,Ns);
             for j=1:Ns
               cosa(j) = sampleHG(ampfunc,dim,1);
             end

        
        
    end

    if sct_type>1
        w=zeros(dim,Ns);
        Nl=size(ow,2);
        q=zeros(1,Ns);
        p=zeros(1,Ns);
        
        for j=1:Ns
           jj=max(ceil(rand*Nl),1);
 
           w(:,j) = rotateBytheta(ow(:,jj),cosa(j),1);
        
           %tp=evalampfunc_general(w(:,j)'*ow,sct_type,ampfunc,dim);
           %tp=tp.^2;
           tp=evaluateHG(w(:,j)'*ow,ampfunc,1,dim);
           
           
           
           p(j)=prod(tp).^(1/Nl);
           q(j)=p(j)./mean(tp);
           %q1(j)=(prod(tp).^(1/Nl));
           %q2(j)=mean(tp);
           
         
        end
        q=q/sum(q);
        cq=cumsum(q);
        r=rand;
        ii=min(find(r<=cq));
        w=w(:,ii);
        pw=p(ii).^0.5;
    end

end