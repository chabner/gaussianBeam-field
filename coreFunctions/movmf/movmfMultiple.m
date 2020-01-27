function [movmf] = movmfMultiple(movmf1,movmf2,matrixMultiple,conjMultiple)
% multiple two movmf components
%
% INPUT
% movmf1,movmf2 : mixture of von Mises-Fisher structres
% matrixMultiple: flag for the result, if matrixMultiple is true, the
% result of movmf is [N1 x N2] vmf mixtures, and if matrixMultiple is
% false, the result is N mixtrs, where N = N1 = N2.
% conjMultiple: active conjegate while multiplying
%
% OUTPUT
% movmf: mixture of vms after multipication, k equals to k1 * k2, and N is
% according to the matrixMultiple flag

if(nargin < 3)
    matrixMultiple = ~(movmf1.N == movmf2.N);
end

if(nargin < 4)
    conjMultiple = true;
end

dim = movmf1.dim;

if(conjMultiple)
    if(~matrixMultiple)
        N = movmf1.N;

        alpha1 = movmf1.alpha;                       % alpha1 in size of [N,k1,1]
        alpha2 = permute(movmf2.alpha,[1,3,2]);      % alpha2 in size of [N,1,k2]
        movmf.alpha = alpha1 .* conj(alpha2);        % movmf.alpha in size of [N,k1,k2]
        movmf.alpha = reshape(movmf.alpha,N,[]);     % movmf.alpha in size of [N,k1 x k2]

        c1 = movmf1.c;                               % c1 in size of [N,k1,1]
        c2 = permute(movmf2.c,[1,3,2]);              % c2 in size of [N,1,k2]
        movmf.c = c1 + conj(c2);                     % movmf.c in size of [N,k1,k2]
        movmf.c = reshape(movmf.c,N,[]);             % movmf.c in size of [N,k1 x k2]

        kappaMu1 = permute(movmf1.mu,[1,2,4,3]);     % kappaMu1 in size of [N,k1,1,dim]
        kappaMu2 = permute(movmf2.mu,[1,4,2,3]);     % kappaMu2 in size of [N,1,k2,dim]
        kappaMu = kappaMu1 + conj(kappaMu2);         % kappaMu in size of [N,k1,k2,dim]
        movmf.mu = reshape(kappaMu,N,[],dim);        % mu in size of [N,k1 x k2,dim]
        

        movmf.k = size(movmf.alpha,2);
        movmf.dim = dim;
        movmf.N = size(movmf.alpha,1);
    else
        N1 = movmf1.N;
        N2 = movmf2.N;

        alpha1 = permute(movmf1.alpha,[1,3,2]);      % alpha1 in size of [N1,1 ,k1,1 ]
        alpha2 = permute(movmf2.alpha,[3,1,4,2]);    % alpha2 in size of [1 ,N2,1 ,k2]
        movmf.alpha = alpha1 .* conj(alpha2);        % movmf.alpha in size of [N1,N2,k1,k2]
        movmf.alpha = reshape(movmf.alpha,N1*N2,[]); % movmf.alpha in size of [N1*N2,k1 x k2]

        c1 = permute(movmf1.c,[1,3,2]);              % alpha1 in size of [N1,1 ,k1,1 ]
        c2 = permute(movmf2.c,[3,1,4,2]);            % alpha2 in size of [1 ,N2,1 ,k2]
        movmf.c = c1 + conj(c2);                     % movmf.alpha in size of [N1,N2,k1,k2]
        movmf.c = reshape(movmf.c,N1*N2,[]);         % movmf.alpha in size of [N1*N2,k1 x k2]

        kappaMu1 = permute(movmf1.mu,[1,5,2,4,3]);   % kappaMu1 in size of [N1,1,k1,1,dim]
        kappaMu2 = permute(movmf2.mu,[5,1,4,2,3]);   % kappaMu2 in size of [1,N2,1,k2,dim]
        kappaMu = kappaMu1 + conj(kappaMu2);         % kappaMu in size of [N1,N2,k1,k2,dim]
        movmf.mu = reshape(kappaMu,N1*N2,[],dim);    % mu in size of [N1 x N2,k1 x k2,dim]

        movmf.k = size(movmf.alpha,2);
        movmf.dim = dim;
        movmf.N = size(movmf.alpha,1);
    end
else
    if(~matrixMultiple)
        N = movmf1.N;

        alpha1 = movmf1.alpha;                       % alpha1 in size of [N,k1,1]
        alpha2 = permute(movmf2.alpha,[1,3,2]);      % alpha2 in size of [N,1,k2]
        movmf.alpha = alpha1 .* alpha2;              % movmf.alpha in size of [N,k1,k2]
        movmf.alpha = reshape(movmf.alpha,N,[]);     % movmf.alpha in size of [N,k1 x k2]

        c1 = movmf1.c;                               % c1 in size of [N,k1,1]
        c2 = permute(movmf2.c,[1,3,2]);              % c2 in size of [N,1,k2]
        movmf.c = c1 + c2;                           % movmf.c in size of [N,k1,k2]
        movmf.c = reshape(movmf.c,N,[]);             % movmf.c in size of [N,k1 x k2]

        kappaMu1 = permute(movmf1.mu,[1,2,4,3]);     % kappaMu1 in size of [N,k1,1,dim]
        kappaMu2 = permute(movmf2.mu,[1,4,2,3]);     % kappaMu2 in size of [N,1,k2,dim]
        kappaMu = kappaMu1 + kappaMu2;               % kappaMu in size of [N,k1,k2,dim]
        movmf.mu = reshape(kappaMu,N,[],dim);        % mu in size of [N,k1 x k2,dim]

        movmf.k = size(movmf.alpha,2);
        movmf.dim = dim;
        movmf.N = size(movmf.alpha,1);
    else
        N1 = movmf1.N;
        N2 = movmf2.N;

        alpha1 = permute(movmf1.alpha,[1,3,2]);      % alpha1 in size of [N1,1 ,k1,1 ]
        alpha2 = permute(movmf2.alpha,[3,1,4,2]);    % alpha2 in size of [1 ,N2,1 ,k2]
        movmf.alpha = alpha1 .* alpha2;              % movmf.alpha in size of [N1,N2,k1,k2]
        movmf.alpha = reshape(movmf.alpha,N1*N2,[]); % movmf.alpha in size of [N1*N2,k1 x k2]

        c1 = permute(movmf1.c,[1,3,2]);              % alpha1 in size of [N1,1 ,k1,1 ]
        c2 = permute(movmf2.c,[3,1,4,2]);            % alpha2 in size of [1 ,N2,1 ,k2]
        movmf.c = c1 + c2;                           % movmf.alpha in size of [N1,N2,k1,k2]
        movmf.c = reshape(movmf.c,N1*N2,[]);         % movmf.alpha in size of [N1*N2,k1 x k2]

        kappaMu1 = permute(movmf1.mu,[1,5,2,4,3]);   % kappaMu1 in size of [N1,1,k1,1,dim]
        kappaMu2 = permute(movmf2.mu,[5,1,4,2,3]);   % kappaMu2 in size of [1,N2,1,k2,dim]
        kappaMu = kappaMu1 + kappaMu2;               % kappaMu in size of [N1,N2,k1,k2,dim]
        movmf.mu = reshape(kappaMu,N1*N2,[],dim);    % mu in size of [N1 x N2,k1 x k2,dim]

        movmf.k = size(movmf.alpha,2);
        movmf.dim = dim;
        movmf.N = size(movmf.alpha,1);
    end
end
end

