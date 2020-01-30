function [conv_vmf]=movmfConv_notSure(apertureVmf,scatteringMovmf)
% convolve aperture with scattering function

% INPUT
% apertureVmf: complex vmf aperture. N = 2, and k = 1.
% scatteringMovmf: vmf mixture of the scattering function, mu only in
% directions of [0,0,1] or [0,0,-1]
%
% OUTPUT
% conv_vmf: convolution result. movmf with N = 2.

% !!!! Need to be checked for backscattering !!!!

% if(apertureVmf.N ~= 2)
%     error('Use this convolution only for two light sources!')
% end

if(apertureVmf.k ~= 1)
    error('Aperture is only one mixture')
end

dim = scatteringMovmf.dim;
k = scatteringMovmf.k;
N = apertureVmf.N;

kappaMu = apertureVmf.mu;                   % kappaMu in size of [N,k,dim]
kappaMu_r = real(kappaMu);
kappa_r = sqrt(sum(kappaMu_r.*(conj(kappaMu_r)),3)); % kappa_r in size of [N,k]
mu_r = kappaMu_r ./ kappa_r;                         % mu_r in size of [N,k,dim]

% aperture_kappa = sqrt(sum(kappaMu.*(conj(kappaMu)),3)); % kappa_r in size of [N,k]
% aperture_mu = kappaMu ./ aperture_kappa;                         % mu_r in size of [N,k,dim]

conv_vmf.mu = zeros(N,k,3,class(kappaMu));
conv_vmf.alpha = zeros(N,k,class(kappaMu));
conv_vmf.c = zeros(N,k,class(kappaMu));
conv_vmf.N = apertureVmf.N;
conv_vmf.k = k;
conv_vmf.dim = dim;
s_conv_kappa = zeros(N,k,class(kappaMu));

if(dim == 3)
    % build rotation matrix that rotates mu_r to the north pole
    % calculate cross and dot products
    % R rotating matrix in size [3,3,N]
    C = cat(3,mu_r(:,:,2),-mu_r(:,:,1),0 * mu_r(:,:,3)); % cross with the north pole
    D = mu_r(:,:,end); % dot with north pole
    
    tFactors = all(C == 0,3);

    Z = zeros(3,3,N,class(C));
    Z(1,3,:) =  C(:,:,2);
    Z(2,3,:) = -C(:,:,1);
    Z(3,1,:) = -C(:,:,2);
    Z(3,2,:) =  C(:,:,1);

    Z_squared = 0 * Z;
    Z_squared(1,1,:) = Z(1,3,:) .* Z(3,1,:);
    Z_squared(1,2,:) = Z(1,3,:) .* Z(3,2,:);
    Z_squared(2,1,:) = Z(2,3,:) .* Z(3,1,:);
    Z_squared(2,2,:) = Z(2,3,:) .* Z(3,2,:);
    Z_squared(3,3,:) = Z(1,3,:) .* Z(3,1,:) + Z(2,3,:) .* Z(3,2,:);

    R = repmat(eye(3),[1,1,N]) + Z + ... 
        Z_squared .* permute((1-D)./(C(:,:,1).^2 + C(:,:,2).^2),[2,3,1]);

    if(any(tFactors))
        R(:,:,tFactors) = repmat(permute(sign(D(tFactors)), [3,4,1,2]),[3,3,1]);
        R(1,2,tFactors) = 0; R(1,3,tFactors) = 0;
        R(2,1,tFactors) = 0; R(2,3,tFactors) = 0;
        R(3,1,tFactors) = 0; R(3,2,tFactors) = 0;
    end
    
    % aperture to gaussian    
    pR = permute(R,[3,1,2]);
    rotatedMu = sum(pR .* apertureVmf.mu,3 );
%     rotatedMu = sum(pR .* aperture_mu,3 );
    rotatedMu = permute(rotatedMu,[1,3,2]);
    
%     A1 = apertureVmf.kappa .* rotatedMu(:,:,3);
    A1 = rotatedMu(:,:,3);
    A1 = [A1,A1,0*A1];
    
%     B1 = apertureVmf.kappa .* permute(rotatedMu(:,:,1:2),[1,3,2]);
%     c1 = apertureVmf.kappa .* rotatedMu(:,:,3) + apertureVmf.c;
    
    B1 = permute(rotatedMu(:,:,1:2),[1,3,2]);
    c1 = rotatedMu(:,:,3) + apertureVmf.c;
    s1 = ones(apertureVmf.N,1,class(C));
    
    for vmfNum = 1:1:k
%         kappa = scatteringMovmf.kappa(1,vmfNum,end);  What about backscattering?!
        kappa = abs(scatteringMovmf.mu(1,vmfNum,end));
        alpha = scatteringMovmf.alpha(vmfNum);
        c = scatteringMovmf.c(vmfNum);

        if (kappa ~= 0)
            % convert to gaussian
            A2 = [kappa, kappa, 0];
            B2 = [0, 0];
            c2 = kappa + c;
            s2 = 1;

            % convolution
            sumA1A2 = A1 + A2;
            detSumA1A2 = sumA1A2(:,1).*sumA1A2(:,2);

            invSumA1A2 = 0 * sumA1A2;
            invSumA1A2(:,1) =  sumA1A2(:,2) ./ detSumA1A2;
            invSumA1A2(:,2) =  sumA1A2(:,1) ./ detSumA1A2;
            invSumA1A2(:,3) = sumA1A2(:,3);

            Bdiff = B1 - B2;

            conv_g_A = 0 * A1;
            conv_g_B = 0 * B1;
            conv_g_c = 0 * c1;
            conv_g_s = 0 * s1;

            conv_g_A(:,1) = A2(:,1) .* invSumA1A2(:,1) .* A1(:,1);
            conv_g_A(:,2) = A2(:,2) .* invSumA1A2(:,2) .* A1(:,2);

            conv_g_B(:,1) = B2(:,1) .* invSumA1A2(:,1) .* A1(:,1) + ...
                            B1(:,1) .* invSumA1A2(:,1) .* A2(:,1);
            conv_g_B(:,2) = B2(:,2) .* invSumA1A2(:,2) .* A1(:,2) + ...
                            B1(:,2) .* invSumA1A2(:,2) .* A2(:,2);

            conv_g_c(:,1) = c1(:,1) + c2(:,1) + ...
                            (Bdiff(:,1).^2) .* invSumA1A2(:,1) / 2 + ...
                            (Bdiff(:,2).^2) .* invSumA1A2(:,2) / 2;

            conv_g_s(:,1) = s1(:,1) .* s2(:,1) .* sqrt(((2*pi).^2) ./ detSumA1A2);
        else
            conv_g_A = 0 * A1;
            conv_g_B = 0 * B1;
            
            detA = A1(:,1).*A1(:,2);
            conv_g_c = c + c1 + log(s1 .* sqrt(((2*pi)^2)./detA)) + ...
                ((B1(:,1).^2).*A1(:,2) + (B1(:,2).^2).*A1(:,1))./(2*detA);
            conv_g_s = 1;
        end
        
        % gaussian back to vmf
        conv_kappa = sqrt(abs(conv_g_A(:,1)).^2 + ...
                          abs(conv_g_B(:,1)).^2 + ... 
                          abs(conv_g_B(:,2)).^2 );
        conv_mu = zeros(N,3,class(C));
        
        conv_mu(:,1) = conv_g_B(:,1) ./ conv_kappa;
        conv_mu(:,2) = conv_g_B(:,2) ./ conv_kappa;
        conv_mu(:,3) = conv_g_A(:,1) ./ conv_kappa;
        
        conv_mu(conv_kappa == 0,1) = 0;
        conv_mu(conv_kappa == 0,2) = 0;
        conv_mu(conv_kappa == 0,3) = 1;

        conv_c = conv_g_c + log(conv_g_s) - conv_g_A(:,1);
        
        % rotate back
        if(isa(R,'gpuArray'))
            conv_mu = pagefun(@mldivide,R,permute(conv_mu,[2,3,1]));
            conv_mu = permute(conv_mu,[3,1,2]);
        else
            for lNum = 1:1:N
                conv_mu(lNum,:) = (R(:,:,lNum) \ conv_mu(lNum,:).').';
            end
        end

        % build the vmf
        conv_vmf.mu(:,vmfNum,:) = conv_mu;
%         conv_vmf.kappa(:,vmfNum) = conv_kappa;
        s_conv_kappa(:,vmfNum) = conv_kappa;
        conv_vmf.alpha(:,vmfNum) = alpha;
        conv_vmf.c(:,vmfNum) = conv_c;
    end   
    
    % normlize to ideal result
    % take absolute
    conv_vmf_abs = movmfAbs(conv_vmf);
    
    % maximal direction is mu_r
    w_max = conv_vmf_abs.mu;
    
    % maximal value in maximal direction
%     log_estimated_conv_max = conv_vmf.kappa .* (sum(conv_vmf.mu .* w_max,3)) + conv_vmf.c;
    log_estimated_conv_max = s_conv_kappa .* (sum(conv_vmf.mu .* w_max,3)) + conv_vmf.c;
    
%     C = sqrt(sum((kappaMu + scatteringMovmf.kappa .* w_max).^2,3));
    C = sqrt(sum((kappaMu + abs(scatteringMovmf.mu(:,:,end)) .* w_max).^2,3));
    
    c = apertureVmf.c + scatteringMovmf.c;
    log_exact_conv_max = c+C + log(2*pi) - log(C);
        
    % fix the maximum where it matters
%     validIdx = abs(estimated_conv_max./exact_conv_max) > 0.01 & abs(estimated_conv_max./exact_conv_max) < 100;
%     conv_vmf.c(validIdx) = conv_vmf.c(validIdx) - log(estimated_conv_max(validIdx)) + log(exact_conv_max(validIdx));
    conv_vmf.c = conv_vmf.c - log_estimated_conv_max + log_exact_conv_max;
    
    conv_vmf.mu = s_conv_kappa .* conv_vmf.mu;
else
    error('2D is not implemented yet')
end

end

