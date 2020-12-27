function samples = vsamp_2(center, kappa, n)
% VSAMP returns a set of points sampled from a vMF
% 
% SAMPLES = VSAMP(CENTER, KAPPA, N)   returns a sample of N points
% from a multi-dimensional von Mises Fisher distribution having
% mean direction given by CENTER and concentration parameter given
% by KAPPA.
% 
% CENTER is a column vector
% This implementation is based on the algorithm VM* of A.T.A. Wood
% A. T. A. Wood. Simulation of the von-Mises Fisher
% Distribution. Comm Stat. Comput. Simul. (1994) v. 23 157-164
%
% NOTE: Current Limitation disallows center(1,1) = 0
% 
% See also BETARND, MLE

% d > 1 of course
d = size(center,1);			% Dimensionality
b = -kappa + sqrt(kappa*kappa + 1);
x0 = (1-b)/(1+b);

samples = zeros(d,n);
activeSamples = 1:1:n;
w = activeSamples;

c = kappa*x0 + 2*log(1-x0*x0);

while (numel(activeSamples) > 0)
    N = numel(activeSamples);
    z = betarnd(1 , 1, 1, N);			% z is a beta rand var
    u = rand(1,N);				% u is unif rand var
    w_act = (1 - (1+b)*z)./(1 - (1-b)*z);
    w(activeSamples) = w_act;
    
    t = kappa*w_act + 2*log(1-x0*w_act) - c;
    activeSamples = activeSamples(t < log(u));
end

v = randn(2,n);			% d-1 dim unit vector from
                                    % unif distr on sphere
v = v ./ sqrt(sum(v.^2,1));
samples(1:2,:) = sqrt(1-w.*w).*v;
samples(3,:) = w;

