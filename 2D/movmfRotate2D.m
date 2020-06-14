function [rotatedMovmf] = movmfRotate2D(movmf,w)
    rotatedMovmf = movmf;
    rotatedMovmf.mu1 = movmf.mu2 .* w(1);
    rotatedMovmf.mu2 = movmf.mu2 .* w(2);
end

