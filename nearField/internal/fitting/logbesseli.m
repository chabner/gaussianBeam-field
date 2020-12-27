function [logb] = logbesseli(dim,x)

if(dim == 3)
    logb = x + log(sqrt(2/pi)/2 ./ sqrt(x));

    x_small_idx = real(x)<5;
    y = x(x_small_idx);

    logb(x_small_idx) = log(sqrt(2/pi) * sinh(y) ./ sqrt(y));
    logb(x == 0) = 0;
end

if(dim == 2)
    if(isa(x,'gpuArray'))
        logb = log(gpuArray(besseli(0,gather(x),1))) + abs(real(x));
    else
        logb = log(besseli(0,x,1)) + abs(real(x));
    end
    
end