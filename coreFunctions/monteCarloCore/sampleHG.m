function costheta=sampleHG(ampfunc, dim,smpNum)
    % sample direction by HG distribution
    g = ampfunc.g;
%     forwardWeight = ampfunc.forwardWeight;
    
%     if(rand > forwardWeight)
%         g = -g;
%     end
    
    if(abs(g) < 0.001)
         costheta = 2*rand(1,1,smpNum)-1;
         return
    end

    if(dim == 3)
%         costheta = (1-g*g)/(1-g+2*g*rand);
%         costheta = (1/(2*g)) * ((1+g*g)-costheta*costheta); %Not calculating theta to save time
        
        costheta = (1-g*g)./(1-g+2*g*rand(1,1,1,1,1,smpNum));
        costheta = (1/(2*g)) .* ((1+g*g)-costheta.^2); %Not calculating theta to save time
    end
    if(dim == 2)
        %for 2D http://www.eugenedeon.com/wp-content/uploads/2016/09/hitchhikers_v0.1.3.pdf
        %chapter 2.5, HG in flatland
        % theta_q = 2*atan( (1-g)/(1+g) * tan( pi/2*(1-2*rand) ) );
        theta_q = 2*atan( (1-g)./(1+g) .* tan( pi/2.*(1-2*rand(1,1,1,1,1,smpNum))) );
        costheta = cos(theta_q);
    end

end