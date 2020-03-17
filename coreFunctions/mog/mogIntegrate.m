function [I] = mogIntegrate(mog)
    I_k = 2 * pi * exp(mog.c + (mog.b1.^2 + mog.b2.^2) ./ (2 * mog.a)) ./ (mog.a);
    I = sum(mog.alpha .* I_k,1);
end

