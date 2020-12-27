#include "standardC.h"

ub32 c_sub2ind(const ub32* sub, ub32 sub_size, const ub32* siz)
{
    ub32 *k, ind = 0, i;
    
    k = (ub32*) malloc(sub_size * sizeof(ub32));
    
    k[0] = 1;
    for(i = 1 ; i < sub_size ; i++)
    {
        k[i] = (siz[i-1]) * k[i - 1];
    }
    
    for(i = 0; i < sub_size ; i++)
    {
        if(siz[i] != 1 && siz[i] != 0)
        {
            ind += sub[i] * k[i];
        }
    }
    
    free(k);
    
    return ind;
}

        
void c_ind2sub(ub32* res, const ub32* siz, ub32 sizDims, ub32 ndx)
{
    ub32 *k, vi, vj, i;
    
    k = (ub32*) malloc(sizDims * sizeof(ub32));
    
    if(sizDims > 2)
    {
        k[0] = siz[0];
        for(i = 1 ; i < sizDims ; i++)
        {
            k[i] = siz[i] * k[i - 1];
        }

        for (i = (sizDims - 1) ; i >= 2 ; i--)
        {
            vi = ndx % k[i-1];
            vj = (ndx - vi)/k[i-1];
            res[i] = vj;
            ndx = vi;
        }
    }

    if(sizDims >= 2)
    {
        vi = ndx % siz[0];
        res[1] = (ndx - vi)/siz[0];
        res[0] = vi;
    }
    else 
    {
        res[0] = ndx;
    }
    
    free(k);
}
