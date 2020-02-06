function [w0,w0p] = vmfFirstDirectionSum_sameBeam(movmf,n)

dim = movmf.dim;

if(dim == 3)
    % normalize
    movmf = movmfAbs(movmf);
    movmf = movmfPick(movmf,n);
    
    % square
    movmf.c = 2 * movmf.c;
    movmf.mu = 2 * movmf.mu;
    movmf.alpha = movmf.alpha.^2;
    
    kappa = sqrt(sum(movmf.mu.^2,3));
    mu = movmfAbsMu(movmf.mu);
    
    log_C = (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(kappa);
    log_C(kappa == 0) = - log(4*pi);
    
    c_1_2 = real(movmf.c);
    alpha_1_2 = abs(movmf.alpha);
    
    alpha = alpha_1_2 .* exp(c_1_2 - log_C);
    
    movmf.alpha = alpha ./ sum(alpha,2);
    movmf.c = log_C;
    
    % sample a direction
    smpNum = datasample(1:1:numel(movmf.alpha),1,'Weights',gather(movmf.alpha));
    w0 = vsamp(permute(gather(mu(1,smpNum,:)),[3,1,2]), gather(kappa(smpNum)), 1);
    
    % calculate probability
    w0p = sqrt(movmfPdf(movmf,w0));
    
    w0 = w0(:);
end

end

