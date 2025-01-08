#ifndef GETMPROF_H
#define GETMPROF_H
#include <stdint.h>
#include "../potential.h"

float getMprof(struct potential **p,int64_t num_p,float *cen,float *mr,float *rpf,int nbins);

#endif 
