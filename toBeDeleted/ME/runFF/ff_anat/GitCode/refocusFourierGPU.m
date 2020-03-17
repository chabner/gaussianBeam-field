function rU=refocus(u,v_max,v_stp,l,box_min,box_max,box_stp,lambda,vsign)
%l is assumed to be a 3x1 ellements vector, the function builds the v grid
%Note: this function assumes an equal resolution in x, y axis
%Also the minimal resolution  is related to the v grid width& spacing 
% To overcome this limitation and allow a potentially faster calculation of sub elements in the greed, need to implement replicas in measurment 
dim=size(box_min,1);
theta_v=[-v_max:v_stp:v_max];

box_stp(1:2)=max(box_stp(1:2)); %set equal resolution

Nv=length(theta_v);
Nv_targ=lambda/2/box_stp(1)/v_stp*2+1;

%theta_v=gpuArray(theta_v);

if Nv_targ<Nv
    error('insufficient target resolution')
end
Nv_padd_size=(Nv_targ-Nv)/2;


gz=[box_min(3):box_stp(3):box_max(3)];
Nz=length(gz);

rU=zeros(Nv_targ,Nv_targ,Nz,'gpuArray');
[vx,vy]=ndgrid(theta_v,theta_v);
vz=sqrt(1-vx.^2-vy.^2)*vsign;
[gx,gy]=ndgrid([box_min(1):box_stp(1):box_max(1)],[box_min(2):box_stp(2):box_max(2)] );    
[Nx,Ny]=size(gx);

Nx_padd_size=(Nv_targ-Nx)/2;

%p_grid=[gx(:)';gy(:)';gz(:)'];
elxy=exp(2*pi*i/lambda*(l(1)*gx+l(2)*gy)); 
%u=conj(u);%.*elxy;

u=reshape(u,Nv,Nv);
u=padarray(u,[Nv_padd_size,Nv_padd_size]);



for j=1:Nz
    evz=exp(-2*pi*i/lambda*gz(j)*vz);
    elz=exp(2*pi*i/lambda*gz(j)*l(3));
    %tu=flipud(fliplr(u));
    tu=u;
    %size(tu)
    %size(tu(1+Nv_padd_size:end-Nv_padd_size,1+Nv_padd_size:end-Nv_padd_size))
    tu(1+Nv_padd_size:end-Nv_padd_size,1+Nv_padd_size:end-Nv_padd_size)=tu(1+Nv_padd_size:end-Nv_padd_size,1+Nv_padd_size:end-Nv_padd_size).*evz;
    rU(:,:,j)=elz*fftshift(fft2(ifftshift(tu)));
    
end
rU=rU/size(rU,1);


rU=rU(Nx_padd_size+1:end-Nx_padd_size,Nx_padd_size+1:end-Nx_padd_size,:);

for j=1:Nz
   rU(:,:,j)=rU(:,:,j).*elxy;
end
return


if (dim==2)
     [gx,gz]=ndgrid([box_min(1):box_stp(1):box_max(1)],[box_min(2):box_stp(2):box_max(2)] );     
     
     p_grid=[gx(:)';gz(:)'];
     
     dimv=size(gx);
     
     vz=sqrt(1-theta_v(:).^2)';
 else
     [gx,gy,gz]=ndgrid([box_min(1):box_stp(1):box_max(1)],[box_min(2):box_stp(2):box_max(2)], [box_min(3):box_stp(3):box_max(3)] );
     
     p_grid=[gx(:)';gy(:)';gz(:)'];
    dimv=size(gx);
    
    [vx,vy]=ndgrid(theta_v,theta_v);
    vz=sqrt(1-vx(:).^2-vy(:).^2)';
    %v=[vx(:)';vy(:)'; vz]; 
    
end
 
Ns=50;

N=size(p_grid,2)

rU=zeros(1,N);

for j=1:Ns:N
    %j/N
    jj=[j+1:min(j+Ns,N)];
    e=exp(-2*pi*i/lambda*(v'*p_grid(:,jj)));
    el=exp(2*pi*i/lambda*(l'*p_grid(:,jj)));
    rU(jj)=(u(:)'*e).*el;
end

rU=reshape(rU,dimv);
