clear

figure
hold on
load('gatheredRes.mat')
A = squeeze(sum(sum(abs(Cs_ff).^2)));
plot(l_base,A/A(1))

load('MEcor_Jan192020_od1.mat')
plot(l(1,:),ncC)