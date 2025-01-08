#ifndef _IO_GADGET_H_
#define _IO_GADGET_H_
#include <stdint.h>
#include "../particle.h"

#define GADGET_HEADER_SIZE 256

// ZXY 2022.02.06: This file is written for loading particle type 2 ......

struct gadgeth2_header {
  uint32_t num_particles[6];
  double   particle_masses[6];
  double   scale_factor;
  double   redshift;
  int32_t  flag_sfr; /* ZXY 2022.02.06: Type 2 particle formation flag */
  int32_t  flag_feedback; /* ZXY 2022.02.06: Type2 particle formation feedback */
  uint32_t num_total_particles[6];
  int32_t  flag_cooling; /* Cooling flag */
  int32_t  num_files_per_snapshot;
  double   box_size;
  double   omega_0;
  double   omega_lambda;
  double   h_0;
  int32_t  flag_stellarage;
  int32_t  flag_metals;
  int32_t  num_total_particles_hw[6]; /* High word of total particle number */
  int32_t  flag_entropy_ics;
  char unused[GADGET_HEADER_SIZE - (sizeof(double)*12) - 
	      (sizeof(int32_t)*13) - (sizeof(uint32_t)*12)];
};  

void load_particles_gadgeth2(char *filename, struct particle **p, int64_t *num_p);

#endif /* _IO_GADGET_H_ */
