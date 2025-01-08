#ifndef _PROPERTIES_H_
#define _PROPERTIES_H_

#define VMAX_BINS 50
#include "io/io_generic.h"

float max_halo_radius_sg(struct halo *h);
void _populate_mass_bins_sg(struct halo *h, struct halo *cur_h, int64_t *bins, int64_t num_bins, float r_scale, int64_t children);
float _estimate_vmax_sg(int64_t *bins, int64_t num_bins, float r_scale);
void estimate_vmax_sg(struct halo *h);
void calc_basic_halo_props_sg(struct halo *h);
void add_ang_mom_sg(double L[3], float c[6], float pos[6]);
void _calc_num_child_particles_sg(struct halo *h);
void calculate_corevel_sg(struct halo *h, struct potential *po, int64_t total_p);
void calc_shape_sg(struct halo *h, int64_t total_p, int64_t bound);
float estimate_total_energy_sg(int64_t total_p, float *energy_ratio);
void _calc_pseudo_evolution_masses_sg(struct halo *h, int64_t total_p, int64_t bound);
void _calc_additional_halo_props_sg(struct halo *h, int64_t total_p, int64_t bound);
void calc_additional_halo_props_sg(struct halo *h);

#endif
