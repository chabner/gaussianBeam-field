function e = evalphaseatt(x,v,sigt,lambda,box_min,box_max,vsign,dirv)
    % complex transmission and attenuation term in multiple scattering
    pv = x'*v;
    bdv = cubeDist(x,box_min,box_max,-vsign * dirv);
    e = exp(1i*pi*2/lambda*pv-sigt/2*bdv);
end