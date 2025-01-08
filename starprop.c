#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <inttypes.h>
#include "rockstar.h"
#include "halo.h"
#include "fof.h"
#include "particle.h"
#include "groupies.h"
#include "subhalo_metric.h"
#include "check_syscalls.h"
#include "config_vars.h"
#include "universal_constants.h"
#include "potential.h"
#include "nfw.h"
#include "distance.h"
#include "fun_times.h"
#include "jacobi.h"
#include "hubble.h"

#define POTENTIAL_ERR_TOL 1.0
#define POTENTIAL_USE_BH 1
#ifndef POTENTIAL_HALT_AFTER_BOUND
#define POTENTIAL_HALT_AFTER_BOUND 0
#endif /* !def POTENTIAL_HALT_AFTER_BOUND */

#define VMAX_BINS 50
#define FAST3TREE_DIM 6
#define POINTS_PER_LEAF 40
#define FAST3TREE_PREFIX SINGLE
#define FAST3TREE_TYPE struct particle
#include "fast3tree.c"

#include "io/io_generic.h"
#include "potential.h"
#include "myutils.h"

int64_t *particle_halos_star = NULL;
float *particle_r_star = NULL;
int64_t num_alloc_pc_star = 0, num_copies_star = 0;

int64_t num_halos_star = 0;
struct halo *halos_star = NULL;
struct extra_halo_info *extra_info_star = NULL;

int64_t min_dens_index_star = 0;
double dynamical_time_star = 0;

int dist_compare_star(const void *a, const void *b) {
  float c = ((struct potential *)a)->r2;
  float d = ((struct potential *)b)->r2;
  if (c>d) return 1;
  if (c<d) return -1;
  return 0;
}

int pot_compare_star(const void *a, const void *b) {
  float c = - ((struct potential *)a)->pe;
  float d = - ((struct potential *)b)->pe;
  if (c>d) return 1;
  if (c<d) return -1;
  return 0;
}

// DY modified, rvir->r200, 18*3.14*3.14~177
double vir_density_star(double a) {
  double x = (Om/pow(a,3))/pow(hubble_scaling(1.0/a-1.0),2.0) - 1.0;
  //return ((18*M_PI*M_PI + 82.0*x - 39*x*x)/(1.0+x));
  return 200.0;
}

//################################################################

double particle_thresh_dens_star[5] = {0};
double particle_rvir_dens_star = 0; //vir_density_star(1.0) * cons_star;
double particle_rvir_dens_star_z0 = 0; // particle_rvir_dens_star;

float _calc_mass_definition_star(char **md) {
  int64_t length = strlen(*md);
  char last_char = (length) ? md[0][length-1] : 0;
  float matter_fraction = (Om/pow(SCALE_NOW,3))/pow(hubble_scaling(1.0/SCALE_NOW-1.0),2.0);
  float cons = Om * CRITICAL_DENSITY / PARTICLE_MASS; // background density
  char *mass = *md;
  float thresh_dens;
  h0=0.671;
  if (mass[0] == 'm' || mass[0] == 'M') mass++;

  if (last_char == 'b' || last_char == 'B')
    thresh_dens = atof(mass) * cons;
  else if (last_char == 'c' || last_char == 'C')
    thresh_dens = atof(mass) * cons / matter_fraction;
  else {
    if (strcasecmp(*md, "vir")) *md = "vir";
    thresh_dens = vir_density_star(SCALE_NOW) * cons;
  }
  particle_rvir_dens_star_z0 = vir_density_star(1.0) * cons;
  //printf(" @@@@@@@@@@ %f\n",thresh_dens/(200.0*CRITICAL_DENSITY/PARTICLE_MASS));
  return thresh_dens;
}

void calc_mass_definition_star(void) {
  //char *vir = "vir";
  char *vir = "200c";
  int64_t i;
  particle_thresh_dens_star[0] = _calc_mass_definition_star(&MASS_DEFINITION); // vir
  particle_thresh_dens_star[1] = _calc_mass_definition_star(&MASS_DEFINITION2); // 200b
  particle_thresh_dens_star[2] = _calc_mass_definition_star(&MASS_DEFINITION3); // 200c
  particle_thresh_dens_star[3] = _calc_mass_definition_star(&MASS_DEFINITION4); // 500c
  particle_thresh_dens_star[4] = _calc_mass_definition_star(&MASS_DEFINITION5); // 2500c
  particle_rvir_dens_star = _calc_mass_definition_star(&vir); // vir is default
  particle_rvir_dens_star_z0 = particle_rvir_dens_star;
  //printf("DY particle_rvir_dens_star = %f",particle_rvir_dens_star);
  dynamical_time_star = 1.0/sqrt((4.0*M_PI*Gc/3.0)*particle_rvir_dens_star*PARTICLE_MASS);
  min_dens_index_star = 0;
  for (i=1; i<5; i++)
    if (particle_thresh_dens_star[i] < particle_thresh_dens_star[min_dens_index_star])
      min_dens_index_star = i;
}


//################################################################

float max_halo_radius_star(struct halo *h) {
  if (LIGHTCONE) lightcone_set_scale(h->pos);
  float thresh_dens = particle_thresh_dens_star[min_dens_index_star]*PARTICLE_MASS;
  float m = (min_dens_index_star) ? h->alt_m[min_dens_index_star-1] : h->m;
  return(cbrt((3.0/(4.0*M_PI))*m/thresh_dens)*1e3);
}


void _populate_mass_bins_star(struct halo *h, struct halo *cur_h, int64_t *bins, int64_t num_bins, float r_scale, int64_t children) {
  int64_t i, j, child, first_child, bin;
  float ds, dx;
  for (i=0; i<cur_h->num_p; i++) {
    for (ds=0,j=0; j<3; j++) {
      dx = h->pos[j]-pot[cur_h->p_start+i].pos[j];
      ds+=dx*dx;
    }
    bin = sqrt(ds)*r_scale;
    if (bin >= num_bins) continue;
    bins[bin]++;
  }

  if (!children) return;
  first_child = child = extra_info_star[cur_h-halos_star].child;
  while (child > -1) {
    _populate_mass_bins_star(h, halos_star + child, bins, num_bins, r_scale, 1);
    child = extra_info_star[child].next_cochild;
    assert(child != first_child);
  }
}

float _estimate_vmax_star(int64_t *bins, int64_t num_bins, float r_scale) {
  int64_t i, tp=0;
  float vmax=0, vcirc, r;
  for (i=0; i<num_bins; i++) {
    r = (i+1.0)/r_scale;
    if (r<FORCE_RES) r=FORCE_RES;
    tp += bins[i];
    vcirc = tp/r;
    if (vcirc > vmax) vmax = vcirc;
  }
  return sqrt(Gc*vmax*PARTICLE_MASS/SCALE_NOW);
}

void estimate_vmax_star(struct halo *h) {
  int64_t bins[VMAX_BINS]={0};
  h->vmax = h->vmax_r = 0;
  if (!(h->child_r>0)) return;
  float r_scale = ((double)VMAX_BINS)/h->child_r;
  _populate_mass_bins_star(h,h,bins,VMAX_BINS,r_scale,0);
  h->vmax = _estimate_vmax_star(bins,VMAX_BINS,r_scale);
  h->vmax_r = sqrt(SCALE_NOW) * _estimate_vmax_star(bins,VMAX_BINS,r_scale)
    * dynamical_time_star;
}

void calc_basic_halo_props_star(struct halo *h) {
  int64_t j, k;
  double pos[6] = {0}, pos2[6] = {0}, x;
  double pos_err, vel_err;
  h->r = h->vrms = 0;
  for (j=0; j<h->num_p; j++)
    for (k=0; k<6; k++) pos[k] += pot[h->p_start + j].pos[k];

  // averaged pos over all halo particles
  for (k=0; k<6; k++) pos[k] /= (double)h->num_p;

  // variance[6]
  for (j=0; j<h->num_p; j++)
    for (k=0; k<6; k++) {
      x = pot[h->p_start + j].pos[k] - pos[k];
      pos2[k] += x*x;
    }

  // sum sigma^2/N
  for (k=0; k<6; k++) {
    if (k<3) h->r += pos2[k] / (double)h->num_p;
    else h->vrms += pos2[k] / (double)h->num_p;
  }

  pos_err = h->r / (double)h->num_p;
  vel_err = h->vrms / (double)h->num_p;

  //if ((!h->min_pos_err) || (h->min_pos_err > pos_err)) {
  //  h->min_pos_err = pos_err;
  //  h->n_core = h->num_p;
  //  for (k=0; k<3; k++) h->pos[k] = pos[k];
  //}

  //if ((!h->min_vel_err) || (h->min_vel_err > vel_err)) {
  //  h->min_vel_err = vel_err;
  //  for (k=3; k<6; k++) h->pos[k] = pos[k];
  //}
  for (k=3; k<6; k++) h->bulkvel[k-3] = pos[k];

  h->m = h->num_p;
  if (!h->num_child_particles) h->num_child_particles = h->num_p;
  
  h->r = cbrt(h->num_p/((4.0*M_PI/3.0)*particle_rvir_dens_star)); // halo rvir
  h->child_r = cbrt(h->num_child_particles/((4.0*M_PI/3.0)*particle_rvir_dens_star));
  //printf("@@@@@@@@@@@@@@@@@@@@@ \nDY DEBUG: %f %f %f\n",h->r,h->num_p,particle_rvir_dens_star);
  estimate_vmax_star(h);
  if (h->vmax_r) h->r = h->vmax_r;
  h->vrms = sqrt(h->vrms);
}

void add_ang_mom_star(double L[3], float c[6], float pos[6]) {
  // L = r x p;
#define cross(a,x,y,s) L[a] s (pos[x]-c[x])*(pos[y+3]-c[y+3])
  cross(0,1,2,+=);
  cross(0,2,1,-=);
  cross(1,2,0,+=);
  cross(1,0,2,-=);
  cross(2,0,1,+=);
  cross(2,1,0,-=);
#undef cross
}

void _calc_num_child_particles_star(struct halo *h) {
  int64_t child, first_child;

  if (h->num_child_particles) return;
  h->num_child_particles = h->num_p;
  first_child = child = extra_info_star[h-halos_star].child;
  while (child > -1) {
    _calc_num_child_particles_star(halos_star + child);
    h->num_child_particles += halos_star[child].num_child_particles;
    child = extra_info_star[child].next_cochild;
    assert(child != first_child);
  }
}

void calc_num_child_particles_star(int64_t h_start) {
  int64_t i;
  for (i=h_start; i<num_halos_star; i++) halos_star[i].num_child_particles = 0;
  for (i=h_start; i<num_halos_star; i++) 
    if (!halos_star[i].num_child_particles) _calc_num_child_particles_star(halos_star + i);
}

void calculate_corevel_star(struct halo *h, struct potential *pot, int64_t total_p) {
  //Assumes pot is already sorted.
  int64_t i, j;
  double vel[3]={0};
  int64_t core_max, rvir_max;
  double var[3]={0}, thisvar, bestvar=0;
  double rvir_thresh = particle_rvir_dens_star*(4.0*M_PI/3.0);

  for (j=total_p-1; j>=0; j--)
    if (j*j > (pot[j].r2*pot[j].r2*pot[j].r2)*(rvir_thresh*rvir_thresh)) break;
  rvir_max = j;
  //printf("@@@@@@@@ rvir_max = %d rvir_thresh = %f\n",rvir_max,rvir_thresh);
  if (rvir_max < 1) return;
  for (j=total_p-1; j>=0; j--)
    if (pot[j].r2*100.0 < pot[rvir_max].r2) break;
  core_max = j;
  
  //printf("@@@@@@@@ core_max = %d\n",core_max);
  if (core_max < 100) core_max = 100;

  for (i=0; i<rvir_max; i++) {
    for (j=0; j<3; j++) {
      double delta = pot[i].pos[j+3] - vel[j];
      vel[j] += delta / ((double)(i+1));
      var[j] += delta * (pot[i].pos[j+3]-vel[j]);
    }
    thisvar = (var[0]+var[1]+var[2]);
    if ((i < 10) || (thisvar < bestvar*(i-3)*i)) {
      if (i > 3) bestvar = thisvar / (double)((i-3)*i);
      else bestvar = 0;
      if (i < core_max) {
	h->n_core = i;
	h->min_vel_err = bestvar;
	for (j=0; j<3; j++) h->corevel[j] = vel[j];
      }
      for (j=0; j<3; j++) h->bulkvel[j] = vel[j];
      h->min_bulkvel_err = bestvar;
    }
  }

  // DY: consider comment the following...
  // for (j=0; j<3; j++) h->pos[j+3] = h->corevel[j];
}


void calc_shape_star(struct halo *h, int64_t total_p, int64_t bound) {

  // ZXY 2022.05.06: Revise to calculate eccentricity within different radius ......
  //                 To do this, add an extra structure in halo *h, named h->starec

  double radi1 = 0.5*1e-3;
  double radi2 = 30*1e-3;
  double radi = radi1;
  int64_t nr;
  for(nr=0; nr<=15; nr++){
 // for(int64_t nr=0; nr<=10; nr++){
  int64_t i,j,k,l,iter=SHAPE_ITERATIONS, analyze_p=0, a,b,c;
  float b_to_a, c_to_a, min_r = FORCE_RES*FORCE_RES;
  double mass_t[3][3], orth[3][3], eig[3]={0},  r=0, dr, dr2, weight=0;
  h->b_to_a = h->c_to_a = 0;
  memset(h->A, 0, sizeof(float)*3);

  // ZXY 2022.05.05: Test which function is used ......
  // printf("\n\n\n ZXY: starprop.c is used !!!!!!!!!!!!??????????????????? \n\n\n ");
  
  radi = radi1 * pow(radi2/radi1, (double)nr/15);  //  0.5 kpc per tick
   
  if (!(h->r>0)) { printf("return1"); return;}
  min_r *= 1e6 / (h->r*h->r);

  for (j=0; j<total_p; j++) {
    // ZXY 2022.05.05: Resume the determination of bounded or not ......
    // ZXY 2022.05.06: Revise to calculate eccentricity within different radius ......
    //                 In addition, pot[j].r2 has been converted from Mpc/h to Mpc in main.c ......
    if (bound && ( (pot[j].pe < pot[j].ke) || (pot[j].r2 > radi*radi ) ) ) continue;
    analyze_p++;
  }
  //printf(">>>>>>>>>> total_p %d analyze_p %d \n",total_p,analyze_p);
  //printf("Halo center: x=%f, y=%f, z=%f \n",h[0].pos[0],h->pos[1],h->pos[2]);
   

  if (analyze_p < 3 || !(h->r>0)) {printf("return2"); return;}
  if (analyze_p < iter) iter = analyze_p;

  // printf("\n\n>>>>>>>>>>>ZXY Check: analyze_p = %d, iter = %d, nr = %d, radi = %f <<<<<<<<<<<\n",analyze_p,iter,nr,radi);

  for (i=0; i<3; i++) {
    memset(orth[i], 0, sizeof(double)*3);
    orth[i][i] = 1;
    eig[i] = (h->r*h->r)*1e-6; // convert to Mpc scale
  }
  // DY DEBUG
  //printf(">>>>>>>>>> h->r %f \n",h->r);
  for (i=0; i<10; i++) {
    for (k=0; k<3; k++) memset(mass_t[k], 0, sizeof(double)*3);
    weight=0;
    for (j=0; j<total_p; j++) {

    // ZXY 2022.04.10: Revised, consider only particles within 25 kpc ......
    // ZXY 2022.05.06: Revise to calculate eccentricity within different radius ......
      if (bound && ( (pot[j].pe < pot[j].ke) || (pot[j].r2 > radi*radi) ) ) continue;
      r=0;
      for (k=0; k<3; k++) {
	for (dr=0,l=0; l<3; l++) {
	  dr += orth[k][l]*(pot[j].pos[l]-h->pos[l]);
	}
	r += dr*dr/eig[k];
      }
      if (r < min_r) r = min_r;
      // DY DEBUG, r>1 !!!!!!!!!!!!!!
      if (!(r>0 && r<=1)) continue;
      //printf(">>>>>>>>>> r %f\n",r);
      double tw = (WEIGHTED_SHAPES) ? 1.0/r : 1.0;
      weight+=tw;
      for (k=0; k<3; k++) {
	dr = pot[j].pos[k]-h->pos[k];
	mass_t[k][k] += dr*dr*tw;
	for (l=0; l<k; l++) {
	  dr2 = pot[j].pos[l]-h->pos[l];
	  mass_t[k][l] += dr2*dr*tw;
	  mass_t[l][k] = mass_t[k][l];
	}
      }
    }

    //printf("weight %f",weight);
    if (!weight) {printf("return3"); return;}
    for (k=0; k<3; k++) for (l=0; l<3; l++) mass_t[k][l] /= (double)weight;
    jacobi_decompose(mass_t, eig, orth);
    a = 0; b = 1; c = 2;
    if (eig[1]>eig[0]) { b=0; a=1; }
    if (eig[2]>eig[b]) { c=b; b=2; }
    if (eig[b]>eig[a]) { int64_t t=a; a=b; b=t; }
    if (!eig[a] || !eig[b] || !eig[c]) {printf("return4");return;}
    b_to_a = sqrt(eig[b]/eig[a]);
    c_to_a = sqrt(eig[c]/eig[a]);
    
    // ZXY 2022.05.06: Warning of !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! return !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((fabs(b_to_a-h->b_to_a) < 0.005*h->b_to_a) &&	(fabs(c_to_a-h->c_to_a) < 0.005*h->c_to_a)) break;
    h->b_to_a = (b_to_a > 0) ? b_to_a : 0;
    h->c_to_a = (c_to_a > 0) ? c_to_a : 0;
    r = sqrt(eig[a]);
    for (k=0; k<3; k++) {
      h->A[k] = 1e3*r*orth[a][k];
      eig[k] *= (h->r*h->r*1e-6)/(r*r);
    }
  // ZXY 2022.05.06: Check the convergence of iteration .....
  // if(i%2==0)
  // printf("\n iteration = %d, b/a = %f, c/a = %f\n", i, b_to_a, c_to_a);

  }  // ket of iteration algorithm of calculating eig[k] and b/a, c/a  ......
  h->starec[0][nr] = radi;
  h->starec[1][nr] = b_to_a;
  h->starec[2][nr] = c_to_a;
  h->starec[3][nr] = c_to_a/b_to_a;
  }  // ket of different radius considered ......
  // printf("\n >>>>>>>>> nr = %d ???????????????????\n",nr);
}

float estimate_total_energy_star(int64_t total_p, float *energy_ratio) {
  int64_t i;
  double phi=0, total_phi = 0, ke = 0, r;
  for (i=total_p-1; i>=0; i--) {
    if (pot[i].pe > pot[i].ke) {
      ke += pot[i].ke;
      r = sqrt(pot[i].r2);
      if (r<FORCE_RES) r = FORCE_RES;
      total_phi += PARTICLE_MASS*i/r + phi;
      phi += PARTICLE_MASS/r;
    }
  }
  total_phi /= 2.0; //U = sum pe/2
  *energy_ratio = 0;
  if (total_phi) *energy_ratio = (ke/total_phi);
  return ((ke - total_phi)*PARTICLE_MASS*Gc/SCALE_NOW);
}

void _calc_pseudo_evolution_masses_star(struct halo *h, int64_t total_p, int64_t bound)
{
  int64_t j, num_part = 0, num_part_pe_d = 0;
  double r, r32, max_pe_b = 0;
 
  //Typical: R_s*4.0; Minimum thresh: R_halo/5.0
  double r_pe_d = h->rs*4.0;
  double r_pe_b = 0;
  if (r_pe_d < h->r/5.0) r_pe_d = h->r/5.0;
  r_pe_d *= 1e-3;
  for (j=0; j<total_p; j++) {
    //if (bound && (pot[j].pe < pot[j].ke)) continue;
    num_part++;
    r = sqrt(pot[j].r2);

    r32 = sqrt(r);
    r32 = r32*r32*r32; //r^(3/2)
    if ((double)(num_part*num_part) / r32 > max_pe_b) {
      max_pe_b = (double)(num_part*num_part) / r32;
      r_pe_b = r;
    }
    
    if (r < r_pe_d) num_part_pe_d = num_part;
  }
  //  if (h->m > 1e13) fprintf(stderr, "%f %f\n", r_pe_b*1e3, h->rs);
  h->m_pe_d = num_part_pe_d * PARTICLE_MASS;
  h->m_pe_b = PARTICLE_MASS*pow(max_pe_b, 2.0/3.0)/
    cbrt(4.0*M_PI*particle_rvir_dens_star_z0/3.0);
}

void _calc_additional_halo_props_star(struct halo *h, int64_t total_p, int64_t bound)
{
  int64_t j, k, part_mdelta=0, num_part=0, np_alt[4] = {0},
    np_vir=0, dens_tot=0, parts_avgd = 0, num_part_half = 0;
  double dens_thresh = particle_thresh_dens_star[0]*(4.0*M_PI/3.0);
  double d1 = particle_thresh_dens_star[1]*(4.0*M_PI/3.0);
  double d2 = particle_thresh_dens_star[2]*(4.0*M_PI/3.0);
  double d3 = particle_thresh_dens_star[3]*(4.0*M_PI/3.0);
  double d4 = particle_thresh_dens_star[4]*(4.0*M_PI/3.0);
  double rvir_thresh = particle_rvir_dens_star*(4.0*M_PI/3.0);
  //printf("d1234 %f %f %f %f %f",d1,d2,d3,d4,rvir_thresh);
  double vmax_conv = PARTICLE_MASS/SCALE_NOW;
  double r, circ_v, vmax=0, rvmax=0, L[3] = {0}, Jh, m=0, ds;
  double lin[3]={3};
  double vrms[3]={0}, xavg[3]={0}, vavg[3]={0}, mdiff;
  double cur_dens, rvir, mvir;
  double djOverdT[3]={0},deOverdT=0;
  float numin=0;

  for (j=0; j<total_p; j++) {
    // No need to distinguish bound and unbound for the stellar component
    // both of them are observable!
    if (bound && (pot[j].pe < pot[j].ke)) continue;
    num_part++;
    r = sqrt(pot[j].r2);
    if (r < FORCE_RES) r = FORCE_RES;
    cur_dens = ((double)num_part/(r*r*r));
    
    if (cur_dens > dens_thresh) {
      part_mdelta = num_part;
      dens_tot = j;
    //}
    //if(r<20.e-3){
    //if(r<12.3e-3){ // 12.4 is r200 computed for Mstar
    //if(r<72.3e-3){
      part_mdelta = num_part;
      dens_tot = j;
    }
    //if(abs(j-num_part)>2) printf("XXXXXXXXXXXXXXXX %d %d",j,num_part);
    //part_mdelta = total_p;
    //dens_tot = total_p;

    if (cur_dens > d1) np_alt[0] = num_part;
    if (cur_dens > d2) np_alt[1] = num_part;
    if (cur_dens > d3) np_alt[2] = num_part;
    if (cur_dens > d4) np_alt[3] = num_part;

    if (cur_dens > rvir_thresh) {
      circ_v = (double)num_part/r;
      np_vir = num_part;
      if (part_mdelta && circ_v > vmax) {
	vmax = circ_v;
	rvmax = r;
      }
    }
  }
  //printf(" XXXXXXXXXXXXXXX tot %d, mdelta %d, dens %d\n",total_p,part_mdelta,dens_tot);
  //total_p=4440;
  //part_mdelta=total_p;
  //dens_tot=total_p;

  h->dens_tot=dens_tot;
  for (j=0; j<dens_tot; j++) {
    if (bound && (pot[j].pe < pot[j].ke)) continue;
    add_ang_mom_star(L, h->pos, pot[j].pos);
    if(pot[j].r2<0.005*0.005){ add_ang_mom_star(lin, h->pos, pot[j].pos); numin++; }
    parts_avgd++;
    if (parts_avgd*2 >= part_mdelta && !num_part_half)
      num_part_half = j;
    for (k=0; k<3; k++) { //Calculate velocity and position averages
      xavg[k] += (pot[j].pos[k]-xavg[k])/((double)parts_avgd);
      mdiff = pot[j].pos[k+3]-vavg[k];
      vavg[k] += mdiff/(double)parts_avgd;
      vrms[k] += mdiff*(pot[j].pos[k+3]-vavg[k]);
    }
    //DY: calculate djOverdT, deOverdT here
    if(bound){
       double dut[6];
       dut[0]=pot[j].pos[0];
       dut[1]=pot[j].pos[1];
       dut[2]=pot[j].pos[2];
       dut[3]=getDUT(0,pot[j].pos);
       dut[4]=getDUT(1,pot[j].pos);
       dut[5]=getDUT(2,pot[j].pos);
       add_dut(djOverdT, h->pos, dut);// take into account a minus sign 
       deOverdT+= - ((pot[j].pos[3]-h->pos[3])*dut[3]+(pot[j].pos[4]-h->pos[4])*dut[4]+(pot[j].pos[5]-h->pos[5])*dut[5]);
       //printf("@@@@@@@@@@@ %f %f %f\n",dut[3],dut[4],dut[5]);
       //printf("@@@@@@@@@@@ %f\n",djOverdT[2]*1e-3);
    }
  }

  m = part_mdelta*PARTICLE_MASS;
  if (!bound) h->m = m;
  else h->mgrav = m;
  for (k=0; k<3; k++) vrms[k] = (parts_avgd>0) ? (vrms[k]/parts_avgd) : 0;
  if ((!bound) == (!BOUND_PROPS)) { //Works even if BOUND_PROPS > 1
    h->Xoff = h->Voff = 0;
    for (k=0; k<3; k++) { 
      ds = xavg[k]-h->pos[k]; h->Xoff += ds*ds;
      ds = vavg[k]-h->pos[k+3]; h->Voff += ds*ds;
    }
    h->alt_m[0] = np_alt[0]*PARTICLE_MASS;
    h->alt_m[1] = np_alt[1]*PARTICLE_MASS;
    h->alt_m[2] = np_alt[2]*PARTICLE_MASS;
    h->alt_m[3] = np_alt[3]*PARTICLE_MASS;
    h->Xoff = sqrt(h->Xoff)*1e3;
    h->Voff = sqrt(h->Voff);
    h->vrms = sqrt(vrms[0] + vrms[1] + vrms[2]); 
    h->vmax = VMAX_CONST*sqrt(vmax*vmax_conv);
    h->rvmax = rvmax*1e3;
    h->halfmass_radius = sqrt(pot[num_part_half].r2)*1e3;
    //h->halfmass_radius = sqrt(pot[2220].r2)*1e3;
    
    h->r = cbrt((3.0/(4.0*M_PI))*np_alt[2]/particle_thresh_dens_star[3])*1e3;
    calc_shape_star(h,np_alt[2],bound);
    h->b_to_a2 = h->b_to_a;
    h->c_to_a2 = h->c_to_a;
    memcpy(h->A2, h->A, sizeof(float)*3);

    h->r = cbrt((3.0/(4.0*M_PI))*part_mdelta/particle_thresh_dens_star[0])*1e3;
    calc_shape_star(h,dens_tot,bound);

    rvir = cbrt((3.0/(4.0*M_PI))*np_vir/particle_rvir_dens_star)*1e3;
    mvir = np_vir*PARTICLE_MASS;
    calc_scale_radius(h, m, h->r, h->vmax, h->rvmax, SCALE_NOW, pot, dens_tot, bound);
    for (j=0; j<3; j++){
       h->J[j] = PARTICLE_MASS*SCALE_NOW*L[j];
       // DY
       h->jin[j] = lin[j]/numin;
       if(bound){
          // [djOverdT] = km/sec Mpc/Gyr
          h->djOverdT[j]=(djOverdT[j]/dens_tot)*3.1536e-3/3.0857;
       }
    }
    // [deOverdT] = km^2/sec^2 / Gyr
    if(bound) h->deOverdT=(deOverdT/dens_tot)*3.1536e-3/3.0857;
    h->energy = estimate_total_energy_star(dens_tot, &(h->kin_to_pot));
    Jh = PARTICLE_MASS*SCALE_NOW*sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
    h->spin = (m>0) ? (Jh * sqrt(fabs(h->energy)) / (Gc*pow(m, 2.5))) : 0;
    h->bullock_spin = (m>0) ? (Jh / (mvir*sqrt(2.0*Gc*mvir*rvir*SCALE_NOW/1e3))) : 0;
    _calc_pseudo_evolution_masses_star(h,total_p,bound);
  }
}

//Assumes center + velocity already calculated.
void calc_additional_halo_props_star(struct halo *h) {
  int64_t j, total_p, k;
  double dens_thresh;

  // default 0
  //if (LIGHTCONE) lightcone_set_scale(h->pos);
  dens_thresh = particle_thresh_dens_star[0]*(4.0*M_PI/3.0);
  if (h->num_p < 1) return;
  //total_p = calc_particle_r_staradii(h, h, h->pos, 0, 0, 0);
  total_p = h->num_p;
  //_reset_potentials(h, h, h->pos, total_p, 0, 0);

  // default 0
  //if (BOUND_OUT_TO_HALO_EDGE) {
  //  qsort(pot, total_p, sizeof(struct potential), dist_compare_star);
  //  for (j=total_p-1; j>=0; j--)
  //    if (j*j / (pot[j].r2*pot[j].r2*pot[j].r2) > dens_thresh*dens_thresh) break;
  //  if (total_p) total_p = j+1;
  //}

  //if (total_p>1) compute_potential(pot, total_p);
  //for (j=0; j<total_p; j++) {
  //  if (pot[j].ke < 0) {
  //    total_p--;
  //    pot[j] = pot[total_p];
  //    j--;
  //  }
  //}

  // sort according to distances
  qsort(pot, total_p, sizeof(struct potential), dist_compare_star);
  // core vel calculated with inner populations
  //calculate_corevel_star(h, pot, total_p);
  //if (extra_info_star[h-halos_star].sub_of > -1)
  //  compute_kinetic_energy(pot, total_p, h->corevel, h->pos);
  //else
  //  compute_kinetic_energy(pot, total_p, h->bulkvel, h->pos);

  _calc_additional_halo_props_star(h, total_p, 0);
  _calc_additional_halo_props_star(h, total_p, 1);

  if (analyze_halo_generic != NULL) analyze_halo_generic(h, pot, total_p);
}
