#ifndef _MYUTILS_H_
#define _MYUTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>


extern double arrMOverR[100000];
extern double getMOR(int64_t i);
extern double getDUT(int dct,float *pos);
extern void setMOR();




//void do_single(int argc, char **argv);

//struct particle *original_p = NULL;
//struct bparticle *bp = NULL;
//int64_t num_bp = 0, num_additional_p = 0;
//struct fast3tree *tree = NULL;
//struct fast3tree_results *rockstar_res = NULL;
//char *skip = NULL;
//struct fof *all_fofs = NULL;
//int64_t num_all_fofs = 0, num_bfofs = 0, num_metafofs = 0;
//int64_t num_fofs_tosend = 0;
//int64_t *fof_order = NULL;
//
//struct particle *copies = NULL; //For storing phase-space FOFs
//int64_t *particle_halos = NULL;
//float *particle_r = NULL;
//struct potential *po = NULL;
//int64_t num_alloc_pc = 0, num_copies = 0;
//
//struct fof *subfofs = NULL;
//int64_t num_subfofs = 0, num_alloced_subfofs = 0;
//
//int64_t num_halos = 0;
//struct halo *halos = NULL;
//struct extra_halo_info *extra_info = NULL;
//
//struct fast3tree_results *res = NULL;
//struct fast3tree *phasetree = NULL;
//
//int64_t num_alloc_gh = 0, num_growing_halos = 0;
//struct halo **growing_halos = NULL;
//
//int64_t *halo_ids = NULL;
//int64_t num_alloced_halo_ids = 0;
//
//double particle_thresh_dens[5] = {0}, particle_rvir_dens = 0,
//  particle_rvir_dens_z0 = 0;
//int64_t min_dens_index = 0;
//double dynamical_time = 0;


#endif /* _MYUTILS_H_ */
