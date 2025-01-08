#ifndef GETSIGMA_H
#define GETSIGMA_H
#include <stdint.h>
#include "../potential.h"

float getsigma(struct potential **p,int64_t num_p,float *cen,float *sigma,float rin,float eps);

#endif 
