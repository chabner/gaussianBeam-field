function w=rotateBytheta(ow,costheta,smpNum)
    % get new vector direction by given theta direction. phi direction is
    % sampled uniformly
    sintheta = real(sqrt(1-real(costheta).^2));
    dim = size(ow,1);
    if(dim == 2)
        sintheta = ((rand(1,1,1,1,1,smpNum)>0.5)*2-1).*sintheta;
        w = cat(1,costheta.*ow(1,1,1,1,1,:) + sintheta.*ow(2,1,1,1,:),...
            -sintheta.*ow(1,1,1,1,1,:) + costheta.*ow(2,1,1,1,:));
    end

    if(dim == 3)
        w = zeros(3,1,1,1,1,smpNum);
        ow1Idx = (abs(ow(3,1,1,1,1,:)) > 0.999999999);
        
        phi = pi*2*rand(1,1,1,1,sum(ow1Idx));
        sinphi = sin(phi);
        cosphi = cos(phi);
        sintheta_ow1 = sintheta(1,1,1,1,1,ow1Idx);
        costheta_ow1 = costheta(1,1,1,1,1,ow1Idx);
        w(1,1,1,1,1,ow1Idx) = sintheta_ow1 .* cosphi;
        w(2,1,1,1,1,ow1Idx) = sintheta_ow1 .* sinphi;
        w(3,1,1,1,1,ow1Idx) = costheta_ow1 .* ow(3,1,1,1,1,ow1Idx) ./ abs(ow(3,1,1,1,1,ow1Idx));
        
        phi = pi*2*rand(1,1,1,1,1,sum(~ow1Idx));
        sinphi = sin(phi);
        cosphi = cos(phi);
        sintheta_now1 = sintheta(1,1,1,1,1,~ow1Idx);
        costheta_now1 = costheta(1,1,1,1,1,~ow1Idx);
        ow1 = ow(1,:,:,:,:,~ow1Idx);
        ow2 = ow(2,:,:,:,:,~ow1Idx);
        ow3 = ow(3,:,:,:,:,~ow1Idx);
        temp = sqrt(1 - ow(3,:,:,:,:,~ow1Idx).^2);
        w(1,1,1,1,1,~ow1Idx) = sintheta_now1.*(ow1.*ow3.*cosphi - ow2.*sinphi)./temp + ow1.*costheta_now1;
        w(2,1,1,1,1,~ow1Idx) = sintheta_now1.*(ow2.*ow3.*cosphi + ow1.*sinphi)./temp + ow2.*costheta_now1;
        w(3,1,1,1,1,~ow1Idx) = -sintheta_now1.*cosphi.*temp+ow3.*costheta_now1;
        
%         if (abs(ow(3)) > 0.999999999)
%             w(1) = sintheta * cosphi;
%             w(2) = sintheta * sinphi;
%             w(3) = costheta * ow(3) / abs(ow(3));
%         else
%             temp = sqrt(1 - ow(3)^2);
%             w(1) = sintheta*(ow(1)*ow(3)*cosphi - ow(2)*sinphi)/temp + ow(1)*costheta;
%             w(2) = sintheta*(ow(2)*ow(3)*cosphi + ow(1)*sinphi)/temp + ow(2)*costheta;
%             w(3) = -sintheta*cosphi*temp+ow(3)*costheta;
%         end

    end
end