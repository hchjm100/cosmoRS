#ifndef GETVPROF_H
#define GETVPROF_H
#include <stdint.h>
#include "../potential.h"

float getsigmaprof(struct potential **p,int64_t num_p,float *cen,float *vpfr,float *vpft,float *vpfp,float *rpf,const int nticks); 

//float getVprof(struct potential **p,int64_t num_p,float *cen,float *rho,float *rpf,int nbins);

#endif 
