function [C]=evalSSCorNF(sigt,albedo,box_min,box_max,l,v,is_ff_l,is_ff_v,maxItr,lambda,doCBS,smpFlg,sct_type,ampfunc,ampfunc0,lsign,vsign)

C=0;
dim=size(box_min,1);
%v=v(:,1);
%l=l(:,1);
for itr=1:maxItr
    x=rand(dim,1).*(box_max-box_min)+box_min;

dv=x-v;
dl=x-l;

n_dl=lsign.*sum(dl.^2,1).^0.5;
n_dv=vsign.*sum(dv.^2,1).^0.5;

ndl=dl./repmat(n_dl,dim,1);
ndv=dv./repmat(n_dv,dim,1);

[pdv] = cubeDist(x,box_min,box_max,-ndv);
[pdl] = cubeDist(x,box_min,box_max,-ndl);

for j=1:size(v,2)
a(j)=evalampfunc_general((-ndv(:,j)'*ndl(:,j)),sct_type,ampfunc,dim);
end
%1./(((2*pi)*n_dl).^0.5)
%1./(((2*pi)*n_dv).^0.5)

C=C+(a(1)*conj(a(2)))*1/abs((2*pi)^((dim-1)*2)*prod(n_dl.*n_dv).^0.5)*exp(-sigt/2*(sum(abs(pdv+pdl))))*exp(2*pi/lambda*(pdv(1)-pdv(2)+pdl(1)-pdl(2)));

end
C=C/maxItr;