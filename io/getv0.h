#ifndef GETV0_H
#define GETV0_H
#include <stdint.h>
#include "../potential.h"

// 1D velocity dispersion at rin+eps
float getv0(struct potential **p,int64_t num_p,float *cen,float rin,float eps,float alpha);

#endif 
