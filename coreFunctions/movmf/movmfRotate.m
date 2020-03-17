function [rotatedMovmf] = movmfRotate(movmf,w)
    rotatedMovmf = movmf;
    rotatedMovmf.mu1 = movmf.mu3 .* w(1);
    rotatedMovmf.mu2 = movmf.mu3 .* w(2);
    rotatedMovmf.mu3 = movmf.mu3 .* w(3);
end

