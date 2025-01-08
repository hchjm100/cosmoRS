#ifndef HALO_H
#define HALO_H

#define HALO_FORMAT_REVISION 2

#include <stdint.h>

struct halo {
  int64_t id,dens_tot;
  float pos[6], corevel[3], bulkvel[3];
  float m, r, child_r, vmax_r, mgrav, vmax, rvmax, rs, klypin_rs, vrms,
    J[3], jin[3], energy, spin, alt_m[4], Xoff, Voff, b_to_a, c_to_a, A[3],
    b_to_a2, c_to_a2, A2[3],
    bullock_spin, kin_to_pot, m_pe_b, m_pe_d, halfmass_radius;
  int64_t num_p, num_child_particles, p_start, desc, flags, n_core;
  float min_pos_err, min_vel_err, min_bulkvel_err;
  float sigma0,csigma,mtme;
  double djOverdT[3],deOverdT;
  double starec[4][20];
};

struct extra_halo_info {
  int64_t child, next_cochild, prev_cochild;
  int64_t sub_of, ph;
  float max_metric;
};

#endif /* HALO_H */
