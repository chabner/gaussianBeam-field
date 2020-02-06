function [movmfRes] = movmfBuild(sctType, ampfunc, dim, k, samples, iterations, useGpu)
% activate optimization code to generate movmf, which is closed to desired
% scattering function
%
% INPUT
% [sctType. ampfunc]: define the desired scattering function.
% sctType = 1: is isotropic function
% sctType = 3: hg scattering function
%   ampfunc.g            : hg g parameter
%   ampfunc.forwardWeight: weight for forward vs backward scattering
% dim: 2/3 dims
% k: number of mixture components
% samples: number of samples of the desired scattering function
% iterations: max iterations of algorithm
%
% OUTPUT
% Mixture of k Von Mises-Fisher scattering functions. k might be tiny then
% the desired k
% movmfRes.k: [1] number of components
% movmfRes.N: [1] number of vmf mixtures (here 1)
% movmfRes.dim: [1] number of vmf dims
% movmfRes.mu: [1,k,dim] central directions of the distribution. this
% algorithm retuns only movmfRes.mu(1,:,end) = {1 , -1}.
% movmfRes.kappa: [1,k] kappa component.
% movmfRes.c: [1,k] normalizing factor.
% movmfRes.alpha: [1,k] wegiht of each component.
%
% The pdf of the distribution is defined as:
% movmf(x) = sum_k{alpha_k * exp(kappa_k *(mu_kT * x) + c)}

if(nargin < 4)
    k = 51;
end

if(nargin < 5)
    if(dim == 2)
        samples = 1e4;
    end
    
    if(dim == 3)
        samples = 1e3;
    end
end

if(nargin < 6)
    iterations = 1e4;
end

if(nargin < 7)
    useGpu = false;
end

if(sctType == 4)
    movmfRes = inf;
    return;
end

if(sctType == 2)
    error('Not implemented yet :\')
end

if(sctType == 1)
    sctType = 3;
    g = 0;
end

if(dim == 2)
    theta = linspace(-pi, pi, samples + 1);
    vectors = [sin(theta(:)),cos(theta(:))];
    normFactor = 1;
end

if(dim == 3)
%     [X,Y,Z] = sphere(samples);
%     vectors = [X(:),Y(:),Z(:)];
%     normFactor = sqrt(X(:).^2 + Y(:).^2);
    theta = linspace(-pi, pi, samples + 1);
    theta = theta(:);
    phi = rand(samples+1,1)*2*pi;
    
    X = sin(theta).*sin(phi);
    Y = sin(theta).*cos(phi);
    Z = cos(theta);
    
    vectors = [X,Y,Z];
    normFactor = sqrt(X(:).^2 + Y(:).^2);
end

if(sctType == 3)
    hgVals = sqrt(evaluateHG(vectors(:,end), ampfunc, 1, dim));
    wtheta = normFactor.*hgVals;
end

movmfRes = wmovmf(vectors,wtheta(:),k,iterations,mean(wtheta),useGpu);

end

