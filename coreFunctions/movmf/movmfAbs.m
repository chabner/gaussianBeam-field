function [movmf] = movmfAbs(movmf)
% Get the abs value of the movmf
%
% INPUT
% movmf: the mixture of vmf we take the absolute value
%
% OUTPUT
% movmf: the movmf after we take the absolute value

movmf.mu1 = real(movmf.mu1);
movmf.mu2 = real(movmf.mu2);
movmf.mu3 = real(movmf.mu3);

movmf.c = real(movmf.c);
movmf.alpha = abs(movmf.alpha);

end

