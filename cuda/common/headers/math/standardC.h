#pragma once
#include "standardTypes.h"
#include<stdlib.h>

ub32 c_sub2ind(const ub32* sub, ub32 sub_size, const ub32* siz);

void c_ind2sub(ub32* res, const ub32* siz, ub32 sizDims, ub32 ndx);
