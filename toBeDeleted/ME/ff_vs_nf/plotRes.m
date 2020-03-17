clear

load('res.mat')

figure
subplot(2,1,1)
hold on

Cs_nf = nf.config_OD_1_ff.Cs_nf ./ (nf.config_OD_1_ff.itersNum + nf.config_OD_1_ff.xRep);
A = squeeze(sum(sum(abs(Cs_nf).^2)));
plot(l_base,A/A(1));

Cs_ff = ff.config_OD_1_ff.Cs_ff;
A = squeeze(sum(sum(abs(Cs_ff).^2)));
plot(l_base,A/A(1));

subplot(2,1,2)
hold on

C_nf = nf.config_OD_1_ff.C_nf ./ (nf.config_OD_1_ff.itersNum + nf.config_OD_1_ff.xRep);
A = squeeze(sum(sum(abs(C_nf).^2)));
plot(l_base,A/A(1));

C_ff = ff.config_OD_1_ff.C_ff;
A = squeeze(sum(sum(abs(C_ff).^2)));
plot(l_base,A/A(1));
