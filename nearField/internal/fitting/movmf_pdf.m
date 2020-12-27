function [vectors,pdf] = movmf_pdf(alpha,mu,kappa,samplesNum,uy)

dim = size(mu,2);
vmfTotalNum = size(mu,1);

if(~(nargin == 5))
    if(nargin == 3)
        samplesNum = 100;
    end

    if(dim == 2)
        theta = linspace(0, 2*pi, samplesNum + 1);
        vectors = [sin(theta(:)),cos(theta(:))];
    else
        [X,Y,Z] = sphere(samplesNum);
        vectors = [X(:),Y(:),Z(:)];
    end
else
    ux = samplesNum;
    X = ux(:);
    Y = uy(:);
    Z = sqrt(1-ux.^2-uy.^2);
    vectors = [X(:),Y(:),Z(:)];
end

logCp  = log(alpha) + (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(dim/2-1,kappa); 

if(dim == 2)
  logCp(kappa == 0) = log(alpha(kappa == 0)) - log(2*pi);
else
  logCp(kappa == 0) = log(alpha(kappa == 0)) - log(4*pi);
end

% Cp = (kappa.^(dim/2 - 1)) ./ ((2*pi).^(dim/2) .* besseli(dim/2 - 1,kappa));
pdf = 0;
for vmfNum = 1:1:vmfTotalNum
    pdf = pdf + ...
        exp(logCp(vmfNum) + kappa(vmfNum) .* (mu(vmfNum,:) * vectors') );
end

end

