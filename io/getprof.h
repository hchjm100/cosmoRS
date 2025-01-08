#ifndef GETPROF_H
#define GETPROF_H
#include <stdint.h>
#include "../potential.h"

float getprof(struct potential **p,int64_t num_p,float *cen,float *rho,float *rpf,int nbins);

#endif 
