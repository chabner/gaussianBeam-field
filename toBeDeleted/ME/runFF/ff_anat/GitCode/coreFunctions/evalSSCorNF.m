function [C]=evalSSCorNF(sigt,albedo,box_min,box_max,l,v,is_ff_l,is_ff_v,maxItr,lambda,doCBS,smpFlg,sct_type,ampfunc,ampfunc0,lsign,vsign)

C=0;
dim=size(box_min,1);
v=v(:,1);
l=l(:,1);
for itr=1:maxItr
    x=rand(dim,1).*(box_max-box_min)+box_min;

dv=x-v;
dl=x-l;

n_dl=lsign.*sum(dl.^2,1).^0.5;
n_dv=vsign.*sum(dv.^2,1).^0.5;

ndl=dl./repmat(n_dl,1,size(l,2));
ndv=dv./repmat(n_dv,1,size(v,2));

[pdv] = cubeDist(x,box_min,box_max,-ndv);
[pdl] = cubeDist(x,box_min,box_max,-ndl);

a=evalampfunc_general((-ndv'*ndl),sct_type,ampfunc,dim);
%1./(((2*pi)*n_dl).^0.5)
%1./(((2*pi)*n_dv).^0.5)

C=C+abs(a).^2*1/abs((2*pi)^2*n_dl.*n_dv)*exp(-sigt*(pdv+pdl));

end
C=C/maxItr;