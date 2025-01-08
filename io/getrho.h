#ifndef GETRHO_H
#define GETRHO_H
#include <stdint.h>
#include "../particle.h"

float getrho(struct particle **p,int64_t num_p,float *cen,float rin);

#endif 
