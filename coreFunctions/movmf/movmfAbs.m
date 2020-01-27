function [movmf] = movmfAbs(movmf)
% Get the abs value of the movmf
%
% INPUT
% movmf: the mixture of vmf we take the absolute value
%
% OUTPUT
% movmf: the movmf after we take the absolute value

movmf.mu = real(movmf.mu);
movmf.c = real(movmf.c);
movmf.alpha = abs(movmf.alpha);

end

