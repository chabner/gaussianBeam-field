function funcVals = evaluateHG(theta, g, dim)

    cosTheta = cos(theta);

    if(dim == 3)
        funcVals = (1 / (4*pi)) * (1 - g * g) ./ ((g * g - 2 * g * cosTheta + 1).^1.5);
    end

    if(dim == 2)
        funcVals = (1 / (2*pi)) * (1 - g * g) ./ (g * g - 2 * g * cosTheta + 1);
    end

end