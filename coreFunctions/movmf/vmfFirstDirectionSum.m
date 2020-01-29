function [w0,w0p] = vmfFirstDirectionSum(movmf,~)

dim = movmf.dim;

if(dim == 3)
    % normalize
    movmf = movmfAbs(movmf);
    
    % square
    movmf.c = 2 * movmf.c;
    movmf.mu = 2 * movmf.mu;
    movmf.alpha = movmf.alpha.^2;
    
    kappa = sqrt(sum(movmf.mu.^2,3));
    mu = movmf.mu ./ kappa;
    mu(~isfinite(mu(:,:,1:2))) = 0;
    mu(~isfinite(mu(:,:,3))) = 1;
    
    log_C = (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(kappa);
    log_C(kappa == 0) = - log(4*pi);
    
    c_1_2 = real(movmf.c);
    alpha_1_2 = abs(movmf.alpha);
    
    alpha = alpha_1_2 .* exp(c_1_2 - log_C);
    
    movmf.alpha = alpha ./ sum(alpha,2);
    movmf.c = log_C;
    
    n = randi(movmf.N);
    movmfn = movmfPick(movmf,n);
    
    % sample a direction
    smpNum = datasample(1:1:numel(movmfn.alpha),1,'Weights',gather(movmfn.alpha));
    w0 = vsamp(permute(gather(mu(1,smpNum,:)),[3,1,2]), gather(kappa(smpNum)), 1);
    
    % calculate probability
    w0p = sqrt(mean(movmfPdf(movmf,w0)));
    
    w0 = w0(:);
end

end

