#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_

#include <stdint.h>
#define POTENTIAL_DONT_CALCULATE_FLAG 1

struct potential {
    // todo
    /* pos[0-2]: position of potential ?
     * pos[3-5]: velocity of potential ?
     * r2: (distance between paticle and center point)^2
     * */
  float pos[6], r2;
  double pe; // potential ?
  float ke; // kinetic energy ?
  int32_t type; // todo: particle type ? definition of different typies
  int32_t flags;

  /*The following fields are not included for the main halo finder. */
#ifdef POTENTIAL_COMPARISON
  float pe2, ke2, pe3, ke3;
  float v,r;
#endif /* POTENTIAL_COMPARISON */
//#ifdef CALC_POTENTIALS
// 2023.05.12 ZXY: Use particle id in potential structure now ......
  int64_t id, hid;
//#endif /* CALC_POTENTIALS */
};

void compute_kinetic_energy(struct potential *po, int64_t num_po, float *vel_cen, float *pos_cen);
void compute_potential(struct potential *po, int64_t num_po);

#endif /* _POTENTIAL_H_ */
