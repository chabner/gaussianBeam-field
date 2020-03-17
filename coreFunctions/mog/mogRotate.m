function [rotatedMog] = mogRotate(mog,w)
%     rotatedMog = movmfToMog(movmfRotate(mogToMovmf(mog),w));
    rotatedMog = mog;
    rotatedMog.b1 = mog.a * w(1) + mog.b1;
    rotatedMog.b2 = mog.a * w(2) + mog.b2;
    rotatedMog.c = mog.c -0.5 * mog.a * (w(1)^2 + w(2)^2) - mog.b1 * w(1) - mog.b2 * w(2);
end

