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
significantDirection = sign(mean(W .* vectors(:,end)));

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
  
  if(mod(iteration,100) == 0)
      finishLoopFlag = 0;
      
      while finishLoopFlag~=2
        finishLoopFlag = 0;
        kappaNum1 = 1;
        while (kappaNum1 <= numel(kappa)) && (finishLoopFlag == 0)
            kappaNum2 = 1;
            while (kappaNum2 <= numel(kappa)) && (finishLoopFlag == 0)
                if(kappaNum1 ~= kappaNum2 && ...
                        abs(kappa(kappaNum1) - kappa(kappaNum2)) < epsilon && ...
                        mu(kappaNum1,end) == mu(kappaNum2,end))
                   kappa(kappaNum2) = [];
                   alpha(kappaNum1) = alpha(kappaNum1) + alpha(kappaNum2);
                   alpha(kappaNum2) = [];
                   mu(kappaNum2,:) = [];
                   
                   finishLoopFlag = 1;
                   k = numel(alpha);
                end
                
                kappaNum2 = kappaNum2 + 1;
            end
            kappaNum1 = kappaNum1 + 1;
        end
        
        if(finishLoopFlag == 0)
            finishLoopFlag = 2;
        end
      end
  end
  
  logNormalize  = log(alpha) + (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(dim,kappa); 
  
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
%   mean(probMat(:) - W(:))
  alpha  = sum(probMat);
  
  normMu   = sqrt(sum((probMat'*vectors).^2,2));
  rbar  = gather(normMu.'./alpha);

  kappa = (rbar*dim - rbar.^3)./(1-rbar.^2);
  
  for newtonIter = 1:1:100
      A_k = besseli(dim/2,kappa,1) ./ besseli(dim/2-1,kappa,1);
      A_k_tag = 1 - A_k.^2 - ((dim - 1) ./ kappa) .* A_k;
      kappa = kappa - (A_k - rbar) ./ A_k_tag;
  end
  
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
mu = gather(mu);
alpha = gather(alpha);
c = (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(dim,kappa);

if(dim == 2)
    c(kappa == 0) = - log(2*pi);
end

if(dim == 3)
    c(kappa == 0) = - log(4*pi);
end

meanVals = 0;
for kNum = 1:1:k
    mu_k = mu(kNum,:);
    if(dim == 3)
        meanVals = meanVals + mean(alpha(kNum) * ...
            exp(kappa(kNum) * mu_k * vectors.' + c(kNum)) .* ...
            sqrt(1-vectors(:,end).^2).');
    end
    
    if(dim == 2)
        meanVals = meanVals + mean(alpha(kNum) * ...
            exp(kappa(kNum) * mu_k * vectors.' + c(kNum)));
    end
end

meanVals = gather(meanVals);

alpha = (hgMean/meanVals) * alpha;
kappaMu = kappa .* mu.';

movmf.mu1 = kappaMu(1,:).';
movmf.mu2 = kappaMu(2,:).';

if(dim == 3)
    movmf.mu3 = kappaMu(3,:).';
end

movmf.alpha = alpha(:);
movmf.c = gather(c(:));
movmf.dim = size(movmf.alpha);
movmf.dim(end+1:5) = 1;

gpuDevice([]);

end