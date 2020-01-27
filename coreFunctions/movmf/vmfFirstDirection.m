function [w0,w0p] = vmfFirstDirection(movmf)

dim = movmf.dim;

if(dim == 3)
    movmf1 = movmfPick(movmf,1);
    movmf2 = movmfPick(movmf,2);
    
    vmfMult = movmfMultiple(movmf1,movmf2,false);
    vmfMult = movmfAbs(vmfMult);
    
    log_C = (dim/2-1)*log(vmfMult.kappa) - (dim/2)*log(2*pi) - logbesseli(dim/2-1,vmfMult.kappa);
    log_C(vmfMult.kappa == 0) = - log(4*pi);
    
    c_1_2 = real(vmfMult.c);
    alpha_1_2 = abs(vmfMult.alpha);
    
    alpha = alpha_1_2 .* exp(c_1_2 - log_C);
    
    % normalize
    vmfMult.alpha = alpha / sum(alpha);
    vmfMult.c = log_C;
    
%     smpNum = randsample(numel(vmfMult.alpha),1,true,vmfMult.alpha);
    smpNum = datasample(1:1:numel(vmfMult.alpha),1,'Weights',gather(vmfMult.alpha));
    w0 = vsamp(gather(permute(vmfMult.mu(1,smpNum,:),[3,1,2])), ...
        (gather(vmfMult.kappa(smpNum))), 1); % not working when mu = [0;0;-1]!!
    
    w0p = sqrt(movmfPdf(vmfMult,w0));
    w0 = w0(:);
end

end

