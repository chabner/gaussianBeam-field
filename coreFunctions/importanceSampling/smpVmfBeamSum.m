function [x,px,missX,n] = smpVmfBeamSum(vmfApperture,smpPreprocess,box_min,box_max,smpNum)
% Sample first scatterer and first scattering direction

% INPUT:
% vmfApperture - vmf apperture
% sigma_t - 1/MFP
% box_min,box_max - [1,dim], the box size
% movmfScatterer - movmf of the scattering phsae function, size of [1,k]

% OUTPUT:
% x - [1,dim] - first scatterer position
% px - fist scatter position probability
% missX - number of trials failed in x

if(nargin == 4)
    smpNum = 1;
end
missX = 0;

%% sample first scatterer
Naperture = prod(vmfApperture.dim(2:end));

% sample (x,y)
[mu_r1, mu_r2, mu_r3] = movmfAbsMu(vmfApperture);
P0_x = -imag(vmfApperture.mu1)/(2*pi);
P0_y = -imag(vmfApperture.mu2)/(2*pi);
P0_z = -imag(vmfApperture.mu3)/(2*pi);

x = zeros(3,1,1,1,1,smpNum);
isOutsideList = 1:1:smpNum;

n = zeros(1,1,1,1,1,smpNum);
probType = zeros(1,1,1,1,1,smpNum);

% weighted random number
Var = 1:1:numel(smpPreprocess.alpha);
Odds = smpPreprocess.alpha; 
Odds = cumsum(Odds./sum(Odds));

while(~isempty(isOutsideList))
    outsideNum = numel(isOutsideList);
    % sample one of the apertures
    apertureNum = randi(Naperture,1,1,1,1,1,outsideNum);
    
    mu_r_n1 = mu_r1(apertureNum); mu_r_n2 = mu_r2(apertureNum); mu_r_n3 = mu_r3(apertureNum);
    P0_n_x = P0_x(apertureNum); P0_n_y = P0_y(apertureNum);  P0_n_z = P0_z(apertureNum);
    
    mu_r_n1 = reshape(mu_r_n1,1,1,1,1,1,outsideNum);
    mu_r_n2 = reshape(mu_r_n2,1,1,1,1,1,outsideNum);
    mu_r_n3 = reshape(mu_r_n3,1,1,1,1,1,outsideNum);
    
    P0_n_x = reshape(P0_n_x,1,1,1,1,1,outsideNum);
    P0_n_y = reshape(P0_n_y,1,1,1,1,1,outsideNum);
    P0_n_z = reshape(P0_n_z,1,1,1,1,1,outsideNum);
    
    mixtureNum = zeros(1,1,1,1,1,outsideNum);
    Pz = zeros(1,1,1,1,1,outsideNum);
    r = zeros(1,1,1,1,1,outsideNum);
    probTypeVec = zeros(1,1,1,1,1,outsideNum);
    
    for iter = 1:1:outsideNum
        mixtureNum(iter) = Var(find(Odds>=rand,1,'first'));
        
%         Pz(iter) = slicesample(smpPreprocess.z0(apertureNum(iter)),1,'logpdf',@(x) ...
%             smpPreprocess.smp_logpdf(mixtureNum(iter),apertureNum(iter),x) );
        
%         Pz(iter) = slicesample_2(smpPreprocess.z0(apertureNum(iter)),1,@(t) ...
%             smpPreprocess.smp_logpdf(mixtureNum(iter),apertureNum(iter),t) );

        if(rand > smpPreprocess.tau)
            uRandNum = rand;

            cdf_a = smpPreprocess.smp_cdf_a(mixtureNum(iter),apertureNum(iter));
            cdf_b = smpPreprocess.smp_cdf_b(mixtureNum(iter),apertureNum(iter));
            cdf_c = smpPreprocess.smp_cdf_c(mixtureNum(iter),apertureNum(iter));
            cdf_fx = @(x) smpPreprocess.smp_cdf_fx(mixtureNum(iter),apertureNum(iter),x);

    %         cdffun = @(x) smpPreprocess.smp_cdf(mixtureNum(iter),apertureNum(iter),x) - uRandNum;
    %         Pz(iter) = fzero(cdffun,smpPreprocess.z0(apertureNum(iter)));

            cdffun = @(x) cdf_a .* imag(cdf_b .* (cdf_fx(x) - cdf_c)) - uRandNum;
            Pz(iter) = fzero(cdffun,smpPreprocess.z0(apertureNum(iter)));
            probTypeVec(iter) = 1;
        else
            Pz(iter) = smpPreprocess.cauchyInvCdf(mixtureNum(iter),apertureNum(iter),rand);
            probTypeVec(iter) = 2;
        end
        
        xy = randn([1,2,1,1,1,1]) .* (smpPreprocess.w0(mixtureNum(iter),apertureNum(iter),Pz(iter)));
        r(iter) = sqrt(xy(1,1,1,1,1,:).^2 + xy(1,2,1,1,1,:).^2);
    end
    
    randDirection = randn(1,3,1,1,1,outsideNum);
    pDir = cross(randDirection,cat(2,mu_r_n1,mu_r_n2,mu_r_n3),2);
    pDir = pDir ./ sqrt(sum(pDir.^2,2));
    
    Px = P0_n_x + Pz .* mu_r_n1 + pDir(1,1,1,1,1,:) .* r;
    Py = P0_n_y + Pz .* mu_r_n2 + pDir(1,2,1,1,1,:) .* r;
    Pz = P0_n_z + Pz .* mu_r_n3 + pDir(1,3,1,1,1,:) .* r;
    
%     Px = P0_n_x + pDir(1,1,1,1,1,:) .* r;
%     Py = P0_n_y + pDir(1,2,1,1,1,:) .* r;
%     Pz = P0_n_z + Pz + pDir(1,3,1,1,1,:) .* r;
    
    isInsideCurrent = (Px >= box_min(1) & Px <= box_max(1)) & ...
                      (Py >= box_min(2) & Py <= box_max(2)) & ...
                      (Pz >= box_min(3) & Pz <= box_max(3));
           
    missX = missX + sum(~isInsideCurrent);
    updateList = isOutsideList(isInsideCurrent);
    isOutsideList = isOutsideList(~isInsideCurrent);
    
    x(1,1,1,1,1,updateList) = Px(isInsideCurrent);
    x(2,1,1,1,1,updateList) = Py(isInsideCurrent);
    x(3,1,1,1,1,updateList) = Pz(isInsideCurrent);
    n(1,1,1,1,1,updateList) = apertureNum(isInsideCurrent);
    probType(1,1,1,1,1,updateList) = probTypeVec(isInsideCurrent);
    
end

projected_z = smpPreprocess.project(x(1,1,1,1,1,:),x(2,1,1,1,1,:),x(3,1,1,1,1,:));
% prob_z = smpPreprocess.eval_pdf(projected_z);

prob_z_att = smpPreprocess.eval_pdf_att(projected_z);
prob_z_cauchy = smpPreprocess.eval_pdf_cauchy(projected_z);
prob_z = (probType == 1) .* prob_z_att + (probType == 2) .* prob_z_cauchy;

projected_vec_x = P0_x + projected_z .* mu_r1;
projected_vec_y = P0_y + projected_z .* mu_r2;
projected_vec_z = P0_z + projected_z .* mu_r3;

D = sqrt( ...
    (x(1,1,1,1,1,:) - projected_vec_x).^2 + ...
    (x(2,1,1,1,1,:) - projected_vec_y).^2 + ...
    (x(3,1,1,1,1,:) - projected_vec_z).^2);

w0_all = smpPreprocess.w0_all(projected_z);

% Rayleigh Distribution with b = w0 parameter, divide by 2*pi*D
prob_xy = (1 ./ w0_all.^2) .* exp(-D.^2 ./ (2 .* w0_all.^2)) ./ (2*pi);

px_full = smpPreprocess.alpha .* prob_z .* prob_xy;

px = sum(mean(reshape(px_full,numel(smpPreprocess.alpha),[],1,1,1,smpNum),2),1);

end
