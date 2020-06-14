function plotFitting(n)
% plot the n-th entry of the fitting table

load('vmfCache','vmfCache')

ampfunc = vmfCache(n).ampfunc;
sctType = vmfCache(n).entry.sctType;
dim = vmfCache(n).entry.dim;

%%
if(sctType == 3)
    theta = linspace(-pi, pi, 3e3 );
    theta = theta(:);

    pf = sqrt(evaluateHG(cos(theta), ampfunc, 1, dim));
end

if(sctType == 2)
    pf = ampfunc.evalAmp;
    if(dim == 3)
        theta = linspace(0, pi, numel(pf) );
    end
    if(dim == 2)
        theta = linspace(0, 2*pi, numel(pf) );
    end
    theta = theta(:);
end

if(dim == 3)
    X = sin(theta);
    Y = sin(theta);
    Z = cos(theta);
    vectors = [X(:),Y(:),Z(:)];
end

if(dim == 2)
    X = sin(theta);
    Z = cos(theta);
    vectors = [X(:),Z(:)];
end

movmfVals = movmfPdf(vmfCache(n).movmf,vectors);

figure
if(sctType == 3 && dim == 3)
    subplot(2,1,1)
end
plot(theta,pf,theta,squeeze(movmfVals));
if(sctType == 3)
    legend(['HG, g = ',num2str(ampfunc.g)],[num2str(vmfCache(n).movmf.dim(1)),' mixtures'])
end

if(sctType == 2)
    legend('Target',[num2str(vmfCache(n).movmf.dim(1)),' mixtures'])
end

% xlim([-pi,pi])
title('full')

%%
if(sctType == 3 && dim == 3)
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
    if(sctType == 3)
        legend(['HG, g = ',num2str(ampfunc.g)],[num2str(vmfCache(n).movmf.dim(1)),' mixtures'])
    end

    if(sctType == 2)
        legend('Target',[num2str(vmfCache(n).movmf.dim(1)),' mixtures'])
    end
    title('zoom')
end

end

