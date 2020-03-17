function [mogMult] = mogMultiple(mog1,mog2,conjMultiple)
% multiple two mog components

if(nargin < 3)
    conjMultiple = true;
end


if(conjMultiple)
    mogMult.alpha = mog1.alpha .* conj(mog2.alpha);
    mogMult.c     = mog1.c      + conj(mog2.c)    ;
    mogMult.b1    = mog1.b1     + conj(mog2.b1)   ;
    mogMult.b2    = mog1.b2     + conj(mog2.b2)   ;
    mogMult.a     = mog1.a      + conj(mog2.a)    ;
else
    mogMult.alpha = mog1.alpha .* mog2.alpha      ;
    mogMult.c     = mog1.c      + mog2.c          ;
    mogMult.b1    = mog1.b1     + mog2.b1         ;
    mogMult.b2    = mog1.b2     + mog2.b2         ;
    mogMult.a     = mog1.a      + mog2.a          ;
end

mogMult.dim = max(mog1.dim,mog2.dim);
	
end