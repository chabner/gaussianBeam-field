f = figure;
subplot(2,3,1)
imagesc(abs(corrMat_up(:,:,11)))
title(['abs up. abs->sum: ',num2str(sum(sum(abs(corrMat_up(:,:,11)))))])
colorbar

subplot(2,3,2)
imagesc(abs(corrMat_right(:,:,11)))
title(['abs right. abs->sum: ',num2str(sum(sum(abs(corrMat_right(:,:,11)))))])
colorbar

subplot(2,3,3)
imagesc(abs(corrMat_axis(:,:,11)))
title(['abs axis. abs->sum: ',num2str(sum(sum(abs(corrMat_axis(:,:,11)))))])
colorbar

subplot(2,3,4)
imagesc(real(corrMat_up(:,:,11)))
title(['real up. sum->abs: ',num2str(abs(sum(sum(corrMat_up(:,:,11)))))])
colorbar

subplot(2,3,5)
imagesc(real(corrMat_right(:,:,11)))
title(['real right. sum->abs: ',num2str(abs(sum(sum(corrMat_right(:,:,11)))))])
colorbar

subplot(2,3,6)
imagesc(real(corrMat_axis(:,:,11)))
title(['real axis. sum->abs: ',num2str(abs(sum(sum(corrMat_axis(:,:,11)))))])
colorbar

f.Position = [0,0,1000,800];

figure
subplot(1,2,1)
dPixel = 30;
angles = deg2rad(0:1:360);
middlePixel = 101;
DDx = round(dPixel * cos(angles)) + middlePixel;
DDy = round(dPixel * sin(angles)) + middlePixel;
A = 0 * angles;

for angNum = 1:1:numel(angles)
    A(angNum) = abs(Cs(DDx(angNum),DDy(angNum),11)).^2;
end

plot(rad2deg(angles),A)
hold on
plot(rad2deg(angles(1)),A(1),'g*','MarkerSize',10)
plot(rad2deg(angles(46)),A(46),'k*','MarkerSize',10)
plot(rad2deg(angles(91)),A(91),'m*','MarkerSize',10)
xlim([0,360])

subplot(1,2,2)
imagesc(abs(Cs(:,:,11)).^2)
hold on
plot(DDy,DDx,'r');

plot(DDy(1),DDx(1),'g*','MarkerSize',10);
plot(DDy(46),DDx(46),'k*','MarkerSize',10);
plot(DDy(91),DDx(91),'m*','MarkerSize',10);

xlim([61,141])
axis equal