function [logb] = logbesseli(x)

logb = x + log(sqrt(2/pi)/2 ./ sqrt(x));

x_small_idx = real(x)<5;
y = x(x_small_idx);

logb(x_small_idx) = log(sqrt(2/pi) * sinh(y) ./ sqrt(y));
logb(x == 0) = 0;
