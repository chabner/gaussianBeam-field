function [movmf] = wmovmf(vectors,W,k,iterations,hgMean,useGpu)
% WMOVMF -- Clusters using mixture of vMF distributions
%
% ---------------------------------------------------------
% CLUST = MOVMF(VECTORS, K)
%    VECTORS matrix of data points, each row is a datapoint
%            these are the vectors to be clustered
%    K       number of clusters desired
%
% CLUST = MOVMF(VECTORS, K, TRUTH)
%    VECTORS the input data points
%    K       the number of clusters
%    TRUTH   true cluster label of each data point
%            each entry is in {1,...,k}
% 
% ---------------------------------------------------------
% Author: Arindam Banerjee
% Minor modifications by: Suvrit Sra 
%
% Copyright lies with the author
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
%
% ---------------------------------------------------------

[D,V] = size(vectors);
dim   = V;

% set initial (and final) directions as forward / backward
mu0 = zeros(1,dim);
mu0(end) = 1;
significantDirection = sign(mean(W .* vectors(:,3)));

if(k == 1)
    mu = significantDirection * mu0;
else
    mu = significantDirection * mu0 .* ones(floor(k/2)+1,1);
    mu = [mu;-significantDirection * mu0 .* ones(k-(floor(k/2)+1),1)];
end

% --------------------------------
% Getting the means from spkmeans
% --------------------------------

if(useGpu)
    mu = gpuArray(mu);
    vectors = gpuArray(vectors);
    W = gpuArray(W);
    arrType = 'gpuArray';
else
    arrType = 'double';
end

%----------------------------------------------
% You can cut the code at this point
%----------------------------------------------

% --------------------------------
% movMF iterations
% --------------------------------

diff      = 100;
epsilon   = 0.0001;
value     = 100;
iteration = 1;

alphaMin = 0.01 / k;

if(k == 1)
    kappaMin = 1;
else
    kappaMin = 1;
end

% initializing kappa, alpha
kappa    = kappaMin*abs(randn(1,k,arrType));
kappa_old = kappa;

if(k > 1)
% if(false)
    kappa(1) = 0;
    dcComponent = true;
else
    dcComponent = false;
end

alpha = ones(1,k) / k;

% W = W / mean(W);

logW = log(W);
meanW = mean(W);

% display('Starting main iterations ...');
while (diff > epsilon  && iteration < iterations)
% while (iteration < 10)
  iteration = iteration + 1;
  
  logNormalize  = log(alpha) + (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(kappa); 
  
  if(dim == 2)
      logNormalize(kappa == 0) = log(alpha(kappa == 0)) - log(2*pi);
  else
      logNormalize(kappa == 0) = log(alpha(kappa == 0)) - log(4*pi);
  end
  
  logProbMat    = (vectors * (mu' .* kappa) + logNormalize); 
  logSum        = log(sum(exp(logProbMat),2)); 
					       
  logProbMat = logProbMat + (logW - logSum);
 
  lpmVals = logProbMat(:);
  
  if(mod(iteration,100) == 0)
      oldvalue = value;
      value = sum(sum(lpmVals(~isinf(lpmVals))));
      diff = oldvalue - value;
      
      if((numel(kappa) == numel(kappa_old)) && (sum(abs(kappa - kappa_old)) < epsilon))
          diff = 0;
      end
      kappa_old = kappa;
  end
  
  % updating component parameters
  probMat = exp(logProbMat);
  alpha  = sum(probMat);
  
  normMu   = sqrt(sum((probMat'*vectors).^2,2));
  rbar  = normMu.'./alpha;
    
  kappa = (rbar*dim - rbar.^3)./(1-rbar.^2);
  alpha = alpha/(D*meanW);
  
  if(dcComponent)
    kappa(1) = 0;
  end
  
  if(alpha(1) < alphaMin)
      dcComponent = false;
  end
  
  legalIdx = alpha > alphaMin;

  alpha = alpha(legalIdx);
  kappa = kappa(legalIdx);
  
  muTmp = mu(:);
  legalIdxMu = repmat(legalIdx,[1,size(mu,2)]);
  mu = reshape(muTmp(legalIdxMu),[],size(mu,2));
  
  if(numel(alpha) < k)
      k = numel(alpha);
      diff = 100;
      value = 1e10;
  end
end


%% Build the movmf
kappa = gather(kappa);
mu = permute(gather(mu),[3,1,2]);

movmf.mu = kappa .* mu;
movmf.c = (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(kappa);

if(dim == 2)
    movmf.c(kappa == 0) = - log(2*pi);
else
    movmf.c(kappa == 0) = - log(4*pi);
end

meanVals = 0;

for kNum = 1:1:k
    mu_k = permute(mu(1,kNum,:),[1,3,2]);
    meanVals = meanVals + mean(alpha(kNum) * ...
        exp(kappa(kNum) * mu_k * vectors.' + movmf.c(kNum)) .* ...
        sqrt(1-vectors(:,3).^2).');
end
movmf.alpha = (hgMean/gather(meanVals)) * gather(alpha);
movmf.k = numel(movmf.alpha);
movmf.N = 1;
movmf.dim = dim;

end