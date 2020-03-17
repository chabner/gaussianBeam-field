% clear
% load('res.mat')
% 
% deltaIdx1 = 1;
% deltaIdx2 = 5;
% 
% s_vmf = vmf.config_OD_3_g098;
% s_gaussian = ff.config_OD_3_g098;

%%

maxVal = max(max(abs([s_vmf.C(:,:,deltaIdx1) , s_vmf.C(:,:,deltaIdx2), s_vmf.Cs(:,:,deltaIdx1), s_vmf.Cs(:,:,deltaIdx1) ])));

f = figure;
f.Position = [0,0,800,1200];
subplot(4,3,1)
imagesc(abs(s_vmf.C(:,:,deltaIdx1)) , [0,maxVal]);
axis equal
title('vmf C - idx1')

subplot(4,3,2)
imagesc(abs(s_vmf.C(:,:,deltaIdx2)) , [0,maxVal]);
axis equal
title('vmf C - idx2')

subplot(4,3,3)
imagesc(abs(s_vmf.C(:,:,deltaIdx2)));
axis equal
title('vmf C - idx2')

subplot(4,3,4)
imagesc(abs(s_vmf.Cs(:,:,deltaIdx1)) , [0,maxVal]);
axis equal
title('vmf Cs - idx1')

subplot(4,3,5)
imagesc(abs(s_vmf.Cs(:,:,deltaIdx2)) , [0,maxVal]);
axis equal
title('vmf Cs - idx2')

subplot(4,3,6)
imagesc(abs(s_vmf.Cs(:,:,deltaIdx2)));
axis equal
title('vmf Cs - idx2')

%%

maxVal = max(max(abs([s_ff.C(:,:,deltaIdx1) , s_ff.C(:,:,deltaIdx2), s_ff.Cs(:,:,deltaIdx1), s_ff.Cs(:,:,deltaIdx1) ])));

subplot(4,3,7)
imagesc(abs(s_ff.C(:,:,deltaIdx1)) , [0,maxVal]);
axis equal
title('ff C - idx1')

subplot(4,3,8)
imagesc(abs(s_ff.C(:,:,deltaIdx2)) , [0,maxVal]);
axis equal
title('ff C - idx2')

subplot(4,3,9)
imagesc(abs(s_ff.C(:,:,deltaIdx2)));
axis equal
title('ff C - idx2')

subplot(4,3,10)
imagesc(abs(s_ff.Cs(:,:,deltaIdx1)) , [0,maxVal]);
axis equal
title('ff Cs - idx1')

subplot(4,3,11)
imagesc(abs(s_ff.Cs(:,:,deltaIdx2)) , [0,maxVal]);
axis equal
title('ff Cs - idx2')

subplot(4,3,12)
imagesc(abs(s_ff.Cs(:,:,deltaIdx2)));
axis equal
title('ff Cs - idx2')