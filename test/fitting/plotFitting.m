function plotFitting(n)
% plot the n-th entry of the fitting table

load('vmfCache','vmfCache')

ampfunc.g = vmfCache(n).g;
ampfunc.forwardWeight = vmfCache(n).forwardWeight;

%%
theta = linspace(-pi, pi, 3e3 );
theta = theta(:);

pf = sqrt(evaluateHG(cos(theta), ampfunc, 1, 3));

X = sin(theta);
Y = sin(theta);
Z = cos(theta);
vectors = [X(:),Y(:),Z(:)];

movmfVals = movmfPdf(vmfCache(n).movmf,vectors);

figure
subplot(2,1,1)
plot(theta,pf,theta,squeeze(movmfVals));
legend(['HG, g = ',num2str(ampfunc.g)],[num2str(vmfCache(n).movmf.dim(1)),' mixtures'])
title('full')

%%
theta = linspace(-pi/20, pi/20, 3e3 );
theta = theta(:);

pf = sqrt(evaluateHG(cos(theta), ampfunc, 1, 3));

X = sin(theta);
Y = sin(theta);
Z = cos(theta);
vectors = [X(:),Y(:),Z(:)];

movmfVals = movmfPdf(vmfCache(n).movmf,vectors);

subplot(2,1,2)
plot(theta,pf,theta,squeeze(movmfVals));
legend(['HG, g = ',num2str(ampfunc.g)],[num2str(vmfCache(n).movmf.dim(1)),' mixtures'])
title('zoom')

end

