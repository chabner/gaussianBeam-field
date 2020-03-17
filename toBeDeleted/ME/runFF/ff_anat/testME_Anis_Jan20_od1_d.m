addpath(genpath('../GitCode'))
lambda=0.5;
lambda=1;
g=0.97;


MFP=100;
sigt=1/MFP;
L=100;  

%L=2730*0.1^3; 
%Lind=6;
albedo=1;
doCBS=0;
smpFlg=2;

v_max=0.16;
v_stp=v_max/20;
l_max=0.15;
l_stp=l_max/4;
stp_r=v_stp/l_stp;
theta_v= (-v_max):v_stp:v_max;
theta_l=[0:40:200,300:100:1000,2000:1000:4000]/(20000+100);

Nl0=length(theta_l);
Nv0=length(theta_v);
Nv=Nv0^2;
Nl=Nl0;
vsign=1;
[vx,vy]=ndgrid(theta_v,theta_v);
v=[vx(:)';vy(:)'; vsign*sqrt(1-vx(:).^2-vy(:).^2)'];
vz=reshape(v(3,:),Nv0,Nv0);


l=[theta_l(:)';zeros(1,Nl);sqrt(1-theta_l(:).^2)'];

xsize=L*10;
box_min=[-xsize;-xsize;-L/2];
box_max=[xsize;xsize;L/2];

maxItr=10^3;
maxItr=10^0;
%maxItr=10^5;
%maxItr=10^2;

corrM=zeros(2,2,length(theta_l));
j1=1; j2=length(theta_l);
%j1=13; j2=24;
%j1=1; j2=12;

randnum=round(rand*1000);

maxItrO=10^4;
maxItrO=1;
corrM=zeros(Nv0,Nv0,Nl);
corrM1=zeros(Nv0,Nv0,Nl);

  c_itr=0;
for itr=1:maxItrO
   
    tic
    for j11=1:Nv0
        for j12=1:Nv0
            for j2=1:Nl
               
          
                
                
                tl=l(:,[1,j2]);
                
                tv=[vx(j11,j12);vy(j11,j12);vz(j11,j12)]*[1,1];
                tv(:,2)=tv(:,1)+(tl(:,2)-tl(:,1));
                tv(3,2)=sqrt(1-sum(tv(1:2,2).^2));
                
                
                   
                

                
                
                [Ms,Mm]=MCcov( sigt,albedo, box_min,box_max, tl, tv,1,1,maxItr,lambda,doCBS,smpFlg,4,g,g);
               
                
                corrM(j11,j12,j2)=Ms(1,2,1,2)+Mm(1,2,1,2);
                corrM1(j11,j12,j2)=Ms(1,2,1,2);

            end
        end
    end
    toc
    c_itr=c_itr+maxItr;
    eval(sprintf('save tmpcor/tmp_corrjan20_od1_%d corrM corrM1 c_itr',randnum))
    
end

d=dir('tmpcor/tmp_corrjan20_od1_*')
tcorrM=0;
tcorrM1=0;
t_c_itr=0;
for j=1:length(d)
 eval(sprintf('load tmpcor/%s',d(j).name))
 tcorrM=tcorrM+corrM;
 tcorrM1=tcorrM1+corrM1;
 t_c_itr=t_c_itr+c_itr;
end

corrM=tcorrM;
corrM1=tcorrM1;


maxC=abs(mean(mean(corrM(:,:,1))))^2;

for j=1:Nl
    figure, imagesc(real(corrM(:,:,j)/max(max(abs(corrM(:,:,1))))*2))
    cC(j)=mean(mean(abs(corrM(:,:,j).^2)));
    tW=abs(corrM(:,:,j).^2);
    evr(j)=sqrt(sum((tW(:)/sum(tW(:))).^2))*0.01*maxC
    %imwrite(imresize(corrM(:,:,j)/max(max(abs(corrM(:,:,1))))*2,4,'nearest'),sprintf('tmp_corr_res/Anis6_%d.jpg',j))
end
evr./cC

figure, plot(l(1,:),cC/cC(1))
ncC=cC/cC(1);

save MEcor_Jan192020_od1 ncC l

return
figure, hold on;
for k=1:length(gL)
ncorr=((abs(corrML(:,k))/abs(corrML(1,k))).^1);
plot(mdeg_theta-mdeg_theta(1),ncorr)
end
legend('g=0','g=0.5','g=0.75','g=0.9')



[meanU]=evalMeanUnifCtr(box_min,box_max,l,v,sigt(1),lambda);
max(abs(meanU))


%figure, plot(mdeg_theta,abs(squeeze(corrM(2,1,:))))
return

load MECBres/CBres_Lind1_theta20_theta480.mat
tcorrM=corrM;
load MECBres/CBres_Lind1_theta500_theta1000.mat
corrM(:,:,1:24)=tcorrM(:,:,1:24);

load MECBres/CBres_Lind3_theta1_theta24.mat
tcorrM=corrM;
load MECBres/CBres_Lind3_theta25_theta200.mat
corrM(:,:,1:24)=tcorrM(:,:,1:24);



load MECBres/CBres_Lind5_theta1_theta12.mat
tcorrM=corrM;
load MECBres/CBres_Lind5_theta13_theta24.mat
corrM(:,:,1:12)=tcorrM(:,:,1:12);
tcorrM=corrM;
load MECBres/CBres_Lind5_theta25_theta200.mat
corrM(:,:,1:24)=tcorrM(:,:,1:24);

k=2*pi/lambda;
tt=((theta_l*k*L)./sinh(theta_l*k*L)).^2;
figure, plot(mdeg_theta,((theta_l*k*L)./sinh(theta_l*k*L)).^2)

