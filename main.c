#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "io/io_gadget.h"
#include "io/getcenter.h"
#include "io/getrho.h"
#include "io/getprof.h"
#include "io/getMprof.h"
#include "io/getsigma.h"
#include "io/getv0.h"
#include "io/io_generic.h"
#include "particle.h"
#include "check_syscalls.h"
#include "myutils.h"
#include "potential.h"
#include "groupies.h"
#include "universal_constants.h"

#include "config_part_mass_ZXY.h"
#include "ZXY_module/calc_pot_MW.h"
// Comoving coordinates
// Remember to set the following two params in config.template.h
// GADGET_LENGTH_CONVERSION = 1

#define RUNSINGLE 1
#define DYDEBUG 1
#define PROFILE 1 
#define RUNSTARS 1
// BH for barnes_hut...
#define POTENTIAL_USE_BH 1

extern struct particle *p   = NULL;
extern struct potential *pot = NULL;
int64_t num_p = 0, num_ph=0, num_ps=0;
int64_t i=0,j=0,k=0;

#include "singleprop.c"
#include "starprop.c"

int main(int argc, char *argv[]){


   printf("~~~~~~~~~ ZXY test: phi(10,20,30) = %f \n", phi_MW(0.10,0.20,0.30));
   const int nfiles=500;
   char* filename[nfiles];
   float cen[6];
   float cenp[6];

   // ZXY 2022.02.12: High oscillation of rho profile, try to enlarge radius step size ......
   // ---------------------------------------------------------------------------------------
   double heps=0.0035;
   FORCE_RES =2e-6;
   const float rmin=heps; // kpc
   const float rmax=20;  // kpc
   int nticks=25;
   // ---------------------------------------------------------------------------------------


   double All_G = 4.30091727e-6;


   float profh[nticks-1];
   float profs[nticks-1];
   float Mprofh[nticks-1];
   float Mprofs[nticks-1];
   float Sprofhr[nticks-1];
   float Sprofht[nticks-1];
   float Sprofhp[nticks-1];
   float Vprofhr[nticks-1];
   float Vprofht[nticks-1];
   float Vprofhp[nticks-1];
   float Vprofs[nticks-1];
   float rpf[nticks];
   float rpf2[nticks];
   float rx=rmin;
   double tor,tom;
   double Hz,dtime;
/*
   memset(rpf, 0, sizeof(float)*(nticks));
   memset(rpf2, 0, sizeof(float)*(nticks));
   memset(profh, 0, sizeof(float)*(nticks-1));
   memset(profs, 0, sizeof(float)*(nticks-1));
   memset(Mprofh, 0, sizeof(float)*(nticks-1));
   memset(Mprofs, 0, sizeof(float)*(nticks-1));
   memset(Sprofhr, 0, sizeof(float)*(nticks-1));
   memset(Sprofht, 0, sizeof(float)*(nticks-1));
   memset(Sprofhp, 0, sizeof(float)*(nticks-1));
   memset(Vprofhr, 0, sizeof(float)*(nticks-1));
   memset(Vprofht, 0, sizeof(float)*(nticks-1));
   memset(Vprofhp, 0, sizeof(float)*(nticks-1));
   memset(Vprofs, 0, sizeof(float)*(nticks-1));
*/
   for(int ni=1;ni<nticks;ni++){ rx=rx*pow(rmax/rmin,1./((float)(nticks-1))); rpf[ni]=rx; }
   printf("\n !!!!!!!!!! rpf[0] = %f ~~~~~~~~~~~~~\n", rpf[0]); 
   //printf("%f %f\n",profh[0],profh[nticks-2]);
   //double scale_m1=0.333;
   double scale_m1=0.;
   cen[0]=0.;cen[1]=0.;cen[2]=0.; // a=3 position
   cen[3]=0; cen[4]=0; cen[5]=0;
   memcpy(cenp,cen,sizeof(float)*6);
   int ifile=0;
   int isnap=atoi(argv[3]);


   // 2023.04.07 ZXY: Time of snapshot ......
   double snapshot2time=0.2;// Gyr per snapshot
   double time_snap = isnap*0.2;


#ifndef RUNSINGLE
   for(ifile=isnap;ifile<nfiles;ifile+=1){
   //for(ifile=isnap;ifile>=0;ifile-=1){
#else
   for(ifile=isnap;ifile<isnap+1;ifile++){
#endif
      filename[ifile] = malloc(sizeof(char)*1000);
      struct particle *ploc = NULL;
      struct particle *pdm = NULL;
      struct particle *ps  = NULL;
      struct potential *potloc = NULL;
      struct potential *poth = NULL;
      struct potential *pots = NULL;
      struct halo *h = NULL;
      struct halo *star = NULL;
      sprintf(filename[ifile],"%s/snapshot_%.3d",argv[1],ifile);
      num_ph=0,num_ps=0;


      // ZXY 2022.02.06: Record type 1 particles as pdm1 and type 2 as pdm2 ......
      int64_t num_ph1, num_ph2;
      struct particle *pdm1 = NULL;
      struct particle *pdm2 = NULL;

      load_particles_gadget2(filename[ifile], &pdm1, &num_ph1);

      PART_MASS_ZXY[1] = (double)PARTICLE_MASS;

      // ZXY 2022.02.06: Add type 2 particles ......
      load_particles_gadgeth2(filename[ifile], &pdm2, &num_ph2);
      PART_MASS_ZXY[2] = (double)PARTICLE_MASS;
      num_ph = num_ph1 + num_ph2;
      printf("\n  >>>>>>>>>>>>> ZXY:  numph1 = %d,  numph2 = %d  \n", num_ph1, num_ph2);
      num_ps=0;
#ifdef RUNSTARS
      load_particles_gadgetStar2(filename[ifile], &ps, &num_ps);
      PART_MASS_ZXY[4] = (double)PARTICLE_MASS;
#endif

     printf("\n >>>>>>>>>>>>>>>> ZXY test: Part Mass 1 = %f, Part Mass 2 = %f,  Part Mass 4 = %f!!!!!!!!!!!!\n", PART_MASS_ZXY[1], PART_MASS_ZXY[2],PART_MASS_ZXY[4]);

    // ZXY note 2022/01/17 : At the end of each load_particles functions, there is a rescale function which convert 
    //                       the length unit from kpc to Mpc. The POSITIONS! are comoving coordinates! while the 
    //                       VELOCITIES are Physical Velocity! Read  io/io_gadget.c -> gadget2_rescale_particles()
    //                       for more details ......

      tor = (double)SCALE_NOW/h0; 
      tom = 1.0/h0;
      printf("!!!!!!!!!!!!!!!!!!! SCALE_NOW = %f!!!!!!!!!!!\n", SCALE_NOW); 
      ploc = (struct particle *) check_realloc(ploc, (num_ph+num_ps)*sizeof(struct particle), "Allocating room for particles.");
      p=ploc;
      //for(i=0;i<num_ph;i++){ p[i].id=i; memcpy(p[i].pos, pdm[i].pos, sizeof(float)*6); }
      //for(j=0;j<num_ps;j++){ p[num_ph+j].id=num_ph+j; memcpy(p[num_ph+j].pos, ps[j].pos, sizeof(float)*6); }

      //for(i=0;i<num_ph1;i++){ p[i].id = i; memcpy(p[i].pos, pdm1[i].pos, sizeof(float)*6); }
      //for(k=0;k<num_ph2;k++){ p[num_ph1+k].id = num_ph1+k; memcpy(p[num_ph1+k].pos, pdm2[k].pos, sizeof(float)*6); }
      //for(j=0;j<num_ps;j++){ p[num_ph+j].id=num_ph+j; memcpy(p[num_ph+j].pos, ps[j].pos, sizeof(float)*6); }
      for(i=0; i<num_ph1; i++)
        { 
         if(i%1000000==0)  printf("ZXY: Check pdm1: p_id = %d , pos = %f,  %f ......\n",pdm1[i].id, pdm1[i].pos[0], pdm1[i].pos[1]);
        }

      // ZXY 2022.02.06 Here we must define room of pdm, then use it, rather than revise it directly !!!!!!!!!!!!!!!
      pdm = (struct particle *) check_realloc(pdm, (num_ph)*sizeof(struct particle), "Allocating room for particles.");

      for(i=0;i<num_ph1;i++){ memcpy(pdm[i].pos, pdm1[i].pos, sizeof(float)*6); }
      for(j=0;j<num_ph2;j++){ memcpy(pdm[num_ph1+j].pos, pdm2[j].pos, sizeof(float)*6); }

      int64_t total_p = num_ph + num_ps;


      // Find the a rough center
#ifdef RUNSINGLE
      //getcenter(&pdm,num_ph,cen,0.018*SCALE_NOW); // unit: Mpc
      getcenter(&ps,num_ps,cen,1.13);
      // getcenter(&p,total_p,cen,1.1);
#else
      //getcenter(&ps,num_ps,cen,0.15);
      getcenter(&pdm,num_ph,cen,0.018*SCALE_NOW); // unit: Mpc
      //getcenter(&pdm,num_ph,cen,2*rmax);
      //getcenter(&p,total_p,cen,0.08);
#endif
   
      //#######################################################################
      // Find potential mininum
      //#######################################################################
      // We calculate potential only once, using the center just found. 
      float dx, dy, dz, r2;
      potloc = (struct potential *) check_realloc(potloc, (total_p)*sizeof(struct potential), "Allocating room for potentials.");
      poth = (struct potential *) check_realloc(poth, (num_ph)*sizeof(struct potential), "Allocating room for potentials.");
      pots = (struct potential *) check_realloc(pots, (num_ps)*sizeof(struct potential), "Allocating room for potentials.");
      memset(potloc, 0, sizeof(struct potential)*(total_p));
      memset(poth, 0, sizeof(struct potential)*(num_ph));
      memset(pots, 0, sizeof(struct potential)*(num_ps));

      pot=potloc;

      //*********************** ZXY 2024.02.08:   CAUTION!!!!!!!!!!!!!!!!!!! This is very dangerous
      //  PARTICLE_MASS = 278.3;

      int nuseful=0; // particles far away from a halo is not useful... (Twice the rvir*a)
      int nuseful_pdm = 0; 
      int nuseful_ps = 0; 

   //-----------------------------------------------------------------------------------------
   // ZXY 2022/01/18: Use pdm. and ps. to record particle positions ......
   // ZXY 2022/01/20: A global revise is, most of "total_p" is replaced by nuseful in following content ......
   // ZXY 2023.07.04: Since pdm[i] do not have id, one need to record id with pdm1 and pdm2 respectively .....


      for(j=0; j<num_ph1; j++) {
        r2 = 0;
        for (k=0; k<3; k++) { 
           dx=pdm1[j].pos[k] - cen[k];
           r2+=dx*dx; 
        }
        if(r2==0) printf("Warning!!!!! r2==0");

        if(1){
           pot[nuseful].type=1;
           pot[nuseful].r2 = r2;
           memcpy(pot[nuseful].pos, pdm1[j].pos, sizeof(float)*6);
           pot[nuseful].id = pdm1[j].id;
           nuseful++;
           nuseful_pdm++;
        }
      }

      for(j=0; j<num_ph2; j++) {
        r2 = 0;
        for (k=0; k<3; k++) {
           dx=pdm2[j].pos[k] - cen[k];
           r2+=dx*dx;
        }
        if(r2==0) printf("Warning!!!!! r2==0");

        if(1){
           pot[nuseful].type=2;
           pot[nuseful].r2 = r2;
           memcpy(pot[nuseful].pos, pdm2[j].pos, sizeof(float)*6);
           pot[nuseful].id = pdm2[j].id;
           nuseful++;
           nuseful_pdm++;
        }
      }
/*
//    2024.01.21  ZXY: For DF2, we need to add star particles for potential calculation !!!!!
      for(j=0; j<num_ps; j++) {
        r2 = 0;
        for (k=0; k<3; k++) {
           dx=ps[j].pos[k] - cen[k];
           r2+=dx*dx;
        }
        if(r2==0) printf("Warning!!!!! r2==0");

        if(1){
           pot[nuseful].type=4;
           pot[nuseful].r2 = r2;
           memcpy(pot[nuseful].pos, ps[j].pos, sizeof(float)*6);
           pot[nuseful].id = ps[j].id;
           nuseful++;
           nuseful_ps++;
        }
      }
*/
      printf("\n>>>>>>ZXY: DM record finished,  nuseful = %d,  nuseful_pdm = %d, nuseful_ps = %d \n",nuseful, nuseful_pdm, nuseful_ps); 


   // ZXY 2022/01/18: Record DM and Star particles respectively ......  
      for(j=0; j<num_ps; j++) {
        r2 = 0;
        for (k=0; k<3; k++) {
           dx=ps[j].pos[k] - cen[k];
           r2+=dx*dx;
        }
        if(r2==0) printf("Warning!!!!! r2==0");
        if(1){
           pot[nuseful].type=4;
           pot[nuseful].r2 = r2;
           memcpy(pot[nuseful].pos, ps[j].pos, sizeof(float)*6);
           pot[nuseful].id = ps[j].id;

           // ZXY 2022/01/26: Record stars independently at first, prepare for calculating star center only ......
           // --------------------------------------------------------
           pots[j].type=4;
           pots[j].r2 = r2;
           memcpy(pots[j].pos, ps[j].pos, sizeof(float)*6);
           // --------------------------------------------------------

           nuseful++;
           nuseful_ps++;
        }
      }

    // -----------------------------------------------------------------------------

      printf(">>>>>>ZXY: Star record finished,  nuseful = %d,  nuseful_pdm = %d, nuseful_ps = %d \n",nuseful, nuseful_pdm, nuseful_ps);


//      free(ploc);
//      ZXY 2022/01/13 : Super strange, if ploc is freed, p is too !!!!!!!!



      // need to sort before computing the potential...
      qsort(pot, nuseful, sizeof(struct potential), dist_compare_sg);

      // ZXY 2022/01/26: calculate stars' potential independently (1) ...... 
      qsort(pots, nuseful_ps, sizeof(struct potential), dist_compare_sg);


      // update and compute the potential mininum
      compute_potential(pot,nuseful);

      // ZXY 2022/01/26: calculate stars' potential independently (2)...... 
      compute_potential(pots,nuseful_ps);


      // Sort by potential energy, prepare for calculating new center (potential minimum) ......
      qsort(pot, nuseful, sizeof(struct potential), pot_compare_sg);

      // ZXY 2022/01/26: Sort stars by their own potential ......
      qsort(pots, nuseful_ps, sizeof(struct potential), pot_compare_sg);

      // ZXY 2022/01/18: See the result of sorting ......
      // for(i=0; i<15; i++) { printf("\n >>>>>>>>> ZXY check: r2= %f, pe= %f",  pot[i].r2 ,  pot[i].pe); }


      // mass center of frac% of particles
      float cx=0,cy=0,cz=0;
      float vx=0,vy=0,vz=0;
      float ninner=0;
      float ninner2=0;

      // ZXY 2022/01/24: Count number of halo and star within the innermose 1000 particles ......
      int innerh = 0, inners = 0;
      int64_t nmfrac=(int64_t) (0.01*nuseful);
      // -------------------------------------------------------------------------

      for(i=0; i<300; i++){
         cx+=pot[i].pos[0]; cy+=pot[i].pos[1]; cz+=pot[i].pos[2];
         vx+=pot[i].pos[3]; vy+=pot[i].pos[4]; vz+=pot[i].pos[5];
         ninner++;
         if(pot[i].type == 1) {innerh++;}
         if(pot[i].type == 4) {inners++;}
      }
#ifdef DYDEBUG
      printf("\n>>>>>> num_ph=%d, num_ps=%d, total_p=%d, nuseful=%d\n", num_ph, num_ps, total_p, nuseful);
      //printf(">>>>>> Find center: fraction of particles outside 10 kpc sphere: %f \n",ninner2/ninner); 
      // printf(">>>>>> Density maximum : cx = %f; cy = %f; cz = %f;\n",cen[0],cen[1],cen[2]);
#endif
      cen[0]=cx/ninner; cen[1]=cy/ninner; cen[2]=cz/ninner;
      cen[3]=vx/ninner; cen[4]=vy/ninner; cen[5]=vz/ninner;
#ifdef DYDEBUG
      printf(">>>>>> ZXY: Potential minimum: cx = %f; cy = %f; cz = %f; with velocity = %f, %f, %f\n",cen[0],cen[1],cen[2],cen[3],cen[4],cen[5]);


      float cendm[6];
      for(j=0; j<6; j++) {cendm[j] = cen[j];}

      printf(">>>>>> !!! ZXY: inner particle: ph = %d, ps = %d\n", innerh, inners);

#endif

      // ZXY 2022/01/24: Recalculate center, using only star particles ......

      float cens[6];
      cx = 0; cy = 0; cz = 0; vx = 0; vy = 0; vz = 0;
      ninner = 0;

      for(i=0; ninner < ( (nuseful_ps<100) ? nuseful_ps:100 ); i++){
         {
            cx+=pots[i].pos[0]; cy+=pots[i].pos[1]; cz+=pots[i].pos[2];
            vx+=pots[i].pos[3]; vy+=pots[i].pos[4]; vz+=pots[i].pos[5];
            ninner++;
         }
      }
      cens[0]=cx/ninner; cens[1]=cy/ninner; cens[2]=cz/ninner;
      cens[3]=vx/ninner; cens[4]=vy/ninner; cens[5]=vz/ninner;
      printf(">>>>>> ZXY: Star center: cx = %f; cy = %f; cz = %f; with velocity = %f, %f, %f\n",cens[0],cens[1],cens[2],cens[3],cens[4],cens[5]);


      // 2024.01.21 ZXY: For DF2, it is star that dominates, especially in later epoch ......
      // 2024.01.31 ZXY: For Crater II, use DM center ..................
      // for(i=0; i<6; i++) {cen[i] = cens[i];}

      float corevel[3]={cen[3],cen[4],cen[5]};
      compute_kinetic_energy(pot, nuseful, corevel, cen);


      // Update the r2 using the newly found center ......
      // ------------------------------------------------
      double dvx, dvy, dvz, dv2;

      for(j=0; j<nuseful; j++) {
        dx = (pot[j].pos[0] - cen[0]);
        dy = (pot[j].pos[1] - cen[1]);
        dz = (pot[j].pos[2] - cen[2]);
        r2 = dx*dx + dy*dy + dz*dz;
        pot[j].r2 = r2;
      }

      qsort(pot, nuseful, sizeof(struct potential), dist_compare_sg);

      // 2023.07.05 ZXY: Recursion process to calculate center more precisely ...
      // ------------------------------------------------------------------------ 
      cx = 0; cy = 0; cz = 0; vx = 0; vy = 0; vz = 0;
      int npcenter = 0;

      for(i=0; (i<nuseful) && (npcenter<200); i++){
        if(pot[i].r2 < 5*5*1e-6 && pot[i].ke < pot[i].pe)
         {
            cx+=pot[i].pos[0]; cy+=pot[i].pos[1]; cz+=pot[i].pos[2];
            vx+=pot[i].pos[3]; vy+=pot[i].pos[4]; vz+=pot[i].pos[5];
            npcenter++;
         }
      }

      cen[0]=cx/npcenter; cen[1]=cy/npcenter; cen[2]=cz/npcenter;
      cen[3]=vx/npcenter; cen[4]=vy/npcenter; cen[5]=vz/npcenter;
      printf("\n >>>>>>!!!! ZXY:  center recalculation: npcenter = %d, cx = %f; cy = %f; cz = %f; with velocity = %f, %f, %f\n", npcenter, cen[0],cen[1],cen[2],cen[3],cen[4],cen[5]);

      for(j=0; j<nuseful; j++) {
        dx = (pot[j].pos[0] - cen[0]);
        dy = (pot[j].pos[1] - cen[1]);
        dz = (pot[j].pos[2] - cen[2]);
        r2 = dx*dx + dy*dy + dz*dz;
        pot[j].r2 = r2;
        if(j%500000==0)
        printf("~~~~~~~~~~ ZXY Check: %d, type = %d, r = %f, potential = %f ...\n", j, pot[j].type, 1000*sqrt(r2), - (All_G * PART_MASS_ZXY[pot[j].type] * pot[j].pe / 1000 ));
      }

      qsort(pot, nuseful, sizeof(struct potential), dist_compare_sg);

      char dmpos[2000];
      char stpos[2000];
      char foldername[2000];

//      printf("\n %s \n", argv[1]);
//      sprintf(foldername,argv[1]);
//      printf("\n %s \n",foldername);
//      printf("\nZXY test..................\n");

      sprintf(stpos,"%s/dataProfile/ST_positions_%.3d.txt", argv[1], isnap);
      sprintf(dmpos,"%s/dataProfile/DM_positions_%.3d.txt", argv[1], isnap);
      printf("\n %s \n",dmpos);

      printf("\n\nBound particle posisions with r < %f kpc recorded in %s \n\n", 8.0, dmpos);

      FILE *fdmp   = fopen(dmpos,"w");
      FILE *fstp   = fopen(stpos,"w");

      double TOTAL_PE = 0;
      double TOTAL_KE = 0;
      double SELF_PE = 0;
      double SELF_KE = 0;
      double EXT_PE = 0;
      double BOUND_PE = 0;
      double BOUND_KE = 0;

      for(j=0; j<nuseful; j++) {
        dx = (pot[j].pos[0] - cen[0]);
        dy = (pot[j].pos[1] - cen[1]);
        dz = (pot[j].pos[2] - cen[2]);
        dvx = pot[j].pos[3] - cen[3];
        dvy = pot[j].pos[4] - cen[4];
        dvz = pot[j].pos[5] - cen[5];
        r2 = dx*dx + dy*dy + dz*dz;
        dv2 = dvx*dvx + dvy*dvy + dvz*dvz;
        pot[j].r2 = r2;

        // todo
        // '/1000': convert to kpc?
        SELF_PE += - ( All_G * PART_MASS_ZXY[pot[j].type] * pot[j].pe / 1000 );
        EXT_PE += PART_MASS_ZXY[pot[j].type] * phi_MW(1000*pot[j].pos[0], 1000*pot[j].pos[1], 1000*pot[j].pos[2]);

        SELF_KE += 0.5 * PART_MASS_ZXY[pot[j].type] * dv2;
        double v2_ZXY = 0;
        // todo
        // v2_ZXY: (absolute velocity of the particle)^2
        // dv2: (relative velocity of the particle based on the center point)^2
        for(k=0; k<3; k++){ v2_ZXY += pot[j].pos[k+3]*pot[j].pos[k+3]; }
        TOTAL_KE += 0.5 * PART_MASS_ZXY[pot[j].type] * v2_ZXY;
        // todo
        // pot[j].ke: self kinetic energy
        // pot[j].pe: self potential
        // only for dark matter: BOUND_PE = SELF_PE, BOUND_KE = SELF_KE
        if(pot[j].ke < pot[j].pe){
          BOUND_PE += - ( All_G * PART_MASS_ZXY[pot[j].type] * pot[j].pe / 1000 );
          BOUND_KE += 0.5 * PART_MASS_ZXY[pot[j].type] * dv2;
        }

	    // r^2 < (10pc)^2
        if( (pot[j].type == 1 || pot[j].type == 2) && r2 < 10*10*1e-12 && pot[j].ke < pot[j].pe )
          {fprintf(fdmp,"%f  %f  %d  %f  %f  %f  %f  %f  %f  \n", 1000*sqrt(r2) , 0.5 * dv2 , pot[j].id, dx, dy, dz, dvx, dvy, dvz);}
        // r^2 < (10pc)^2
        else if( (pot[j].type == 4) && r2 < 10*10*1e-6 && pot[j].ke < pot[j].pe )
          {fprintf(fstp,"%f  %f  %d  %f  %f  %f  %f  %f  %f  \n", 1000*sqrt(r2) , 0.5 * dv2 , pot[j].id, dx, dy, dz, dvx, dvy, dvz);}
      }

      TOTAL_PE = SELF_PE + EXT_PE;

      fclose(fdmp);
      fclose(fstp);
      // After this process, r2 should be the distance to the potential minimum ......
      // ------------------------------------------------


      // ZXY 2023.07.04: Count r < 1.066 bound particles to calculate Vc(R_1/2) of CraterII ...
      double npdm_Rh = 0;
      double npdm_10 = 0, npdm_20 = 0, npdm_50 = 0, npdm_100 = 0,  npdm_200 = 0, npdm_500 = 0;
      // 2024.07.13 Calculate toatl potential and kinetic energy ......

      for (j=0; j<nuseful; j++) {
 
        if(pot[j].r2 < 1.066*1.066*1e-6 && pot[j].ke < pot[j].pe && pot[j].type != 4)
	{
          npdm_Rh += 1;
          // todo
          /* distance < 10 pc
           * npdm_10: # of '< 10 pc'
           * */
          if(pot[j].r2 < 0.01*0.01*1e-6){npdm_10 += 1.0;}
	  if(pot[j].r2 < 0.02*0.02*1e-6){npdm_20 += 1.0;}
          if(pot[j].r2 < 0.05*0.05*1e-6){npdm_50 += 1.0;}
          if(pot[j].r2 < 0.1*0.1*1e-6){npdm_100 += 1.0;}
          if(pot[j].r2 < 0.2*0.2*1e-6){npdm_200 += 1.0;}
        }
      }

      double mdm_Rh = PART_MASS_ZXY[1] * npdm_Rh;
      double mdm_10 = PART_MASS_ZXY[1] * npdm_10;
      double mdm_20 = PART_MASS_ZXY[1] * npdm_20;
      double mdm_50 = PART_MASS_ZXY[1] * npdm_50;
      double mdm_100 = PART_MASS_ZXY[1] * npdm_100;
      double mdm_200 = PART_MASS_ZXY[1] * npdm_200;

      printf("\n !!!!!! !!!!!  Inner 1.066 kpc particle number: %f !!!!!!!, mass: %f\n", npdm_Rh, mdm_Rh);
      printf("\n !!!!!! !!!!!  Check particle mass: %f !!!!!!!\n", PARTICLE_MASS);


      float rhoh = 0;
      float rhos = 0;

      memcpy(cenp,cen,sizeof(float)*6);
      rhoh = getrho(&pdm,num_ph,cenp,heps*SCALE_NOW); // compute density in r<0.3 kpc

      printf("\n >>>>>>>> ZXY test: cenpdm = %f, %f, %f, %f, %f, %f ......\n",cenp[0], cenp[1], cenp[2], cenp[3], cenp[4], cenp[5]);


      memcpy(cenp,cen,sizeof(float)*6);
      rhos = getrho(&ps,num_ps,cenp,heps*SCALE_NOW);

      printf(">>>>>>> ZXY test: cenps = %f, %f, %f, %f, %f, %f ...... \n",cenp[0], cenp[1], cenp[2], cenp[3], cenp[4], cenp[5]);
/*
      int ibin = 0;
      double ZXYSx[nticks-1];
      double ZXYSy[nticks-1];
      double ZXYSz[nticks-1];
      double countbin[nticks-1];
      for(ibin=0; ibin<nticks-1; ibin++){
        ZXYSx[ibin] = 0; ZXYSy[ibin] = 0; ZXYSz[ibin] = 0; countbin[ibin] = 0;
      }
      printf("center check: cen = %f, %f, %f, %f, %f, %f~~~~~~~~~~~\n", cen[0], cen[1], cen[2], cen[3], cen[4], cen[5]);
      for(i=0; i<total_p; i++){
	if(pot[i].ke < pot[i].pe && pot[i].type != 4){
          vx = pot[i].pos[3] - cen[3]; 
          vy = pot[i].pos[4] - cen[4]; 
          vz = pot[i].pos[5] - cen[5];
     
          ibin = 0;
          while(sqrt(1e6*pot[i].r2) > rpf[ibin+1] && ibin<(nticks-1)) ibin++;
          if(ibin < nticks-1){
            countbin[ibin]+=1.0;
            ZXYSx[ibin]+=(vx*vx);
            ZXYSy[ibin]+=(vy*vy);
            ZXYSz[ibin]+=vz*vz;
	  }
        }
      }
      for(ibin=0; ibin<nticks-1; ibin++){
        if(countbin[ibin]>=20.0){
          ZXYSx[ibin] = sqrt(ZXYSx[ibin]/countbin[ibin]);
          ZXYSy[ibin] = sqrt(ZXYSy[ibin]/countbin[ibin]);
          ZXYSz[ibin] = sqrt(ZXYSz[ibin]/countbin[ibin]);	
          printf("~~~~   ZXY test DM dispersion: count = %f, sigx0 = %f, sigy0 = %f, sigz0 = %f   ~~~~~~~~~~~~~\n",countbin[ibin] , ZXYSx[ibin],ZXYSy[ibin],ZXYSz[ibin]);
	}
        else{
          ZXYSx[ibin] = 0;
          ZXYSy[ibin] = 0;
          ZXYSz[ibin] = 0;
        }
      }
*/


#ifdef DYDEBUG
      //printf("N halos: %d, N stars: %d\n",num_ph,num_ps);
      printf(">>>>>> Central core density, rhoh = %f, rhos = %f\n",rhoh,rhos);
#endif
       
      j=0;k=0;
      for(i=0; i < nuseful; i++) {
         if(pot[i].type==1){ memcpy(&poth[j],&pot[i],sizeof(struct potential)*1); j++;}
         if(pot[i].type==4){ memcpy(&pots[k],&pot[i],sizeof(struct potential)*1); k++;}
      }
      free(potloc);


      h    = check_realloc(h,sizeof(struct halo),"Halo");
      memset(h,0,sizeof(struct halo));
      memcpy(h->pos, cen, sizeof(float)*6);

      h->corevel[0]=cen[3]; h->corevel[1]=cen[4]; h->corevel[2]=cen[5];

      h->num_p=num_ph; 
      total_p=num_ph;
      pot=poth;
    
      // Remove unbinding particles -------------------------------------
      float cenb[3]={0,0,0};
      float cenunb[3]={0,0,0};
      int64_t nb=0,nunb=0;
      float nh76=0,ns76=0;
      float nh70=0,ns70=0;
   
      const nfreq=1;


      // ZXY 2022/01/21: Resume this block to count bounded and unbounded mass ......
      for (j=0; j<total_p; j++) {
        if ( (pot[j].ke-pot[j].pe) > 0 ) { // unbound particles
             nunb++; 
             cenunb[0]+=pot[j].pos[0];
             cenunb[1]+=pot[j].pos[1];
             cenunb[2]+=pot[j].pos[2];
          total_p--;
          pot[j] = pot[total_p];
          j--;
        }
        else{ // bounded particles
           nb++;
           cenb[0]+=pot[j].pos[0];
           cenb[1]+=pot[j].pos[1];
           cenb[2]+=pot[j].pos[2];
           r2 = 0;
           for(k=0; k<3; k++) { dx = pot[j].pos[k]-cen[k]; r2 += dx*dx; }
        }
      }

/*
      if(ifile == 0) 
      {
         boundmass = fopen("../dataProfile/boundmass.txt","w");
         fprintf(boundmass, "FORCE_RES = %f\n", FORCE_RES);
      }
      else
      {  boundmass = fopen("../dataProfile/boundmass.txt","a+");  }
*/

/*
      // ZXY 2022.01.29: Bound mass should not include all "bounded particles" whose pe > ke, there should be a
      //                 range that exclude far away particles ...... calc_additional_halo_props will set h->m
      //                 or h->mgrav to this value ......
      char bdmass[1000];

      FILE *boundmass;

      sprintf(bdmass,"%s/boundmass.txt",argv[1]);
      boundmass = fopen(bdmass,"a+");

      r2 = sqrt( cen[0]*cen[0] + cen[1]*cen[1] + cen[2]*cen[2] );

      fprintf(boundmass, "%f %f %f %f %f %f\n", SCALE_NOW, PARTICLE_MASS*total_p/1e10, cen[0], cen[1], cen[2], r2);
*/


   #ifdef PROFILE
   //if((ifile+4)%10==0){ // make sure that the last snaphot is included
/*
   int ibin = 0;
   float ZXYSx[nticks-1];
   float ZXYSy[nticks-1];
   float ZXYSz[nticks-1];
   float countbin[nticks-1];
   printf("center check: cen = %f, %f, %f, %f, %f, %f", cen[0], cen[1], cen[2], cen[3], cen[4], cen[5]);
   for(i=0; i<total_p; i++){
     if(abs(pot[i].pos[3])>0){ vx = pot[i].pos[3] - cen[3]; }
     else{
       printf("!!!!!!!!!!!!!!! vx error, potvx = %f, cen[3] = %f !!!!!!!!!!!!!!!!!\n", pot[i].pos[3], cen[3]);
       continue;
     }
     if(abs(pot[i].pos[4])>0){ vy = pot[i].pos[4] - cen[4]; }
     else{
        printf("!!!!!!!!!!!!!!! vy error !!!!!!!!!!!!!!!!!\n");
        continue;
     }
     if(abs(pot[i].pos[5])>0){ vz = pot[i].pos[5] - cen[5]; }
     else{
       printf("!!!!!!!!!!!!!!! vz error !!!!!!!!!!!!!!!!!\n");
       continue;
     }

     ibin = 0;
     while(1e6*pot[i].r2 > rpf[ibin+1] && ibin<(nticks-1)) ibin++;

     countbin[ibin]+=1.0;
     ZXYSx[ibin]+=(vx*vx);
     ZXYSy[ibin]+=(vy*vy);
     ZXYSz[ibin]+=vz*vz;
     if(i%10000000==0)
       printf("~~~~~~~check SY: ZXYSy[0] = %f  ~~~~~~~~~~~~~~\n", ZXYSy[0]);
   }
   for(ibin=0; ibin<nticks-1; ibin++){
     if(countbin[ibin]>20){
       ZXYSx[ibin] = sqrt(ZXYSx[ibin]/countbin[ibin]);
       ZXYSy[ibin] = sqrt(ZXYSy[ibin]/countbin[ibin]);
       ZXYSz[ibin] = sqrt(ZXYSz[ibin]/countbin[ibin]);
       printf("~~~~   ZXY test DM dispersion: count = %f, sigx0 = %f, sigy0 = %f, sigz0 = %f   ~~~~~~~~~~~~~\n",countbin[ibin] , ZXYSx[ibin],ZXYSy[ibin],ZXYSz[ibin]);
     }
     else{
       ZXYSx[ibin] = 0;
       ZXYSy[ibin] = 0;
       ZXYSz[ibin] = 0;
     }
   } 
   printf("~~~~~~~~~~~   ZXY test DM dispersion: sigx0 = %f, sigy0 = %f, sigz0 = %f   ~~~~~~~~~~~~~~~~~ \n", ZXYSx[0],ZXYSy[0],ZXYSz[0]);
 */
   

   memcpy(cenp,cen,sizeof(float)*6);

   if(1){
      //Density profile: Mass unit: Solar mass; Length unit: kpc

      getprof(&pot, total_p,cenp,profh,rpf,nticks); // pot includes only the bounded DM
      getMprof(&pot, total_p,cenp,Mprofh,rpf,nticks); // pot includes only the bounded DM

      getsigmaprof(&pot, total_p,cenp,Sprofhr,Sprofht,Sprofhp,rpf,nticks); // pot includes only the bounded DM
      getVprof(&pot, total_p,cenp,Vprofhr,Vprofht,Vprofhp,rpf,nticks); // pot includes only the bounded DM

      char outpfh[1000];
      sprintf(outpfh,"%s/dataProfile/DMprof_snap%.3d.txt", argv[1], ifile);
      FILE *fbh   = fopen(outpfh,"w");
      // ZXY 2022/01/21: If use "a+" as opening mode, contents will be added to the end of the file, which 
      //                 can be used to record halo's bounded mass  ......

      fprintf(fbh,"# r    rho   mr   sigmar  sigmat  sigmap vr vt vp beta\n");
      for(int i=0;i<nticks-1;i++) {
         float VT2=Vprofht[i]*Vprofht[i]+Vprofhp[i]*Vprofhp[i];
         float beta=1.0-VT2/(2*Vprofhr[i]*Vprofhr[i]) ;
         //printf(" rx2 %d %f \n",xi,(float)(rpf2[xi]+rpf2[xi+1])/2.);
         //printf(" rx3 %d %f \n",xi,(float)(rpf[xi]+rpf[xi+1])/2.);
         fprintf(fbh,  "%f    %f    %f   %f   %f   %f   %f   %f   %f   %f  \n",(rpf[i]+rpf[i+1])/2.,profh[i],Mprofh[i],Sprofhr[i],Sprofht[i],Sprofhp[i],Vprofhr[i],Vprofht[i],Vprofhp[i],beta);
      }
      fclose(fbh);


      #ifdef RUNSTARS
      if(num_ps > 0){ // 2024.08.14  ZXY: if no star involved, do not run it to avoid crush ...... (1/3)

      float Vprofsr[nticks-1];
      float Vprofst[nticks-1];
      float Vprofsp[nticks-1];
      float Sprofsr[nticks-1];
      float Sprofst[nticks-1];
      float Sprofsp[nticks-1];

      memcpy(rpf2,rpf,sizeof(float)*(nticks));

//    ZXY 2022.04.20 Set star profile coordinates here !!!!!!!!!!!
/*      float rmax2, rmin2;
      rmax2 = 30;
      rmin2 = 0.6;
      rx = rmin2;
      for(int ni=1;ni<nticks;ni++){ rx=rx*pow(rmax2/rmin2,1./((float)(nticks-1))); rpf2[ni]=rx; }

      getprof(&pots,num_ps ,cenp,profs,rpf2,nticks);
      getMprof(&pots,num_ps ,cenp,Mprofs,rpf2,nticks);
      getsigmaprof(&pots,num_ps,cenp,Sprofsr,Sprofst,Sprofsp,rpf2,nticks);
      getVprof(&pots,num_ps ,cenp,Vprofsr,Vprofst,Vprofsp,rpf2,nticks);


      char outpfs[1000];
      sprintf(outpfs,"%s/dataProfile/prof_star_%s_snap%.3d.txt", argv[1] ,argv[2],ifile);

      FILE *fstar = fopen(outpfs,"w");
      fprintf(fstar,"# r    rho    mr   sigmar   sigmat   sigmap   vr   vt   vp   beta\n");
      for(int i=0;i<nticks-1;i++) {
      fprintf(fstar,"%f    %f    %f    %f    %f    %f    %f    %f    %f\n",(rpf2[i]+rpf2[i+1])/2.,profs[i],Mprofs[i],Sprofsr[i],Sprofst[i],Sprofsp[i],Vprofsr[i],Vprofst[i],Vprofsp[i]);
      }
      fclose(fstar);
*/
      } // End of num_ps > 0 ...... (1/3) 
      #endif
   }
   #endif

      //printf("Nbounded: %d\n",total_p);
      qsort(pot, total_p, sizeof(struct potential), dist_compare_sg);
      h->num_p=total_p; // update halo particle number


#ifdef DYDEBUG
/*
      h->pos[0]=pot[indexmin].pos[0];
      h->pos[1]=pot[indexmin].pos[1];
      h->pos[2]=pot[indexmin].pos[2];
      cen[0]=pot[indexmin].pos[0];
      cen[1]=pot[indexmin].pos[1];
      cen[2]=pot[indexmin].pos[2];
*/   
      printf("------------------ HALO ANALYSIS --------------------\n");
      printf("-------------------------------------------------------\n");
      printf("Number of binding particles: %d, total_p (bounded) = %d .\n",nb,total_p);
      printf("Binding particle fraction:               %f\n",(float)nb/((float)num_ph));
      printf("Halo particle fraction:                  %f\n",(float)total_p/((float)num_ph));
      printf(">>>> Mass center of binding particles:   %f %f %f\n",cenb[0]/((float)total_p),cenb[1]/((float)total_p),cenb[2]/((float)total_p));
      if((num_ph-total_p)>0) printf(">>>> Mass center of unbinding particles: %f %f %f\n",cenunb[0]/((float)(num_ph-total_p)),cenunb[1]/((float)(num_ph-total_p)),cenunb[2]/((float)(num_ph-total_p)));
#endif
   
      calc_mass_definition_sg();
      calc_basic_halo_props_sg(h);
      calc_additional_halo_props_sg(h);
 
      h->sigma0=getv0(&pot,total_p,cen,0.0,heps*SCALE_NOW,0);
      h->csigma=0;

      float sigmaH[3]={0,0,0};
      getsigma(&pot,total_p,cen,sigmaH,0.,heps*SCALE_NOW);


#ifdef DYDEBUG
      //printf("pe[100]: %f\n",pot[100].pe);
      printf("-------------------------------------------------------\n");
      printf("\n----RockStar results for a single halo-----------------\n");
      printf("pos[6]: %f %f %f %f %f %f\n",h->pos[0],h->pos[1],h->pos[2],h->pos[3],h->pos[4],h->pos[5]);
      printf("core vel %f %f %f \n",h->corevel[0],h->corevel[1],h->corevel[2]);
      printf("bulk vel %f %f %f \n",h->bulkvel[0],h->bulkvel[1],h->bulkvel[2]);
      printf("mass: %e\n",h->m); // vir
      printf("mbound: %e\n",h->mgrav);
      printf("r: %f\n",h->r); // vir
      printf("rs: %f\n",h->rs);
      printf("num_p: %d\n",h->num_p);
      printf("dens_tot: %d\n",h->dens_tot);
      printf("vrms: %f\n",h->vrms);
      printf("vmax: %f\n",h->vmax);
      printf("rvmax: %f\n",h->rvmax);
      printf("J[3] %e %e %e\n",h->J[0],h->J[1],h->J[2]);
      printf("j[3] %f %f %f\n",h->J[0]/h->m,h->J[1]/h->m,h->J[2]/h->m);
      printf("jin[3] (r<5 kpc) %f %f %f\n",h->jin[0],h->jin[1],h->jin[2]);
      printf("djOverdT[3] %f %f %f\n",h->djOverdT[0],h->djOverdT[1],h->djOverdT[2]);
      printf("deOverdT: %f\n",h->deOverdT);
      printf("b_to_a: %f\n",h->b_to_a);
      printf("c_to_a: %f\n",h->c_to_a);
      printf("A[3]: %f %f %f\n",h->A[0],h->A[1],h->A[2]);
      printf("kin_to_pot: %f\n",h->kin_to_pot);
      printf("energy: %e\n",h->energy);
      printf("spin: %f\n",h->spin);
      printf("Xoff: %f\n",h->Xoff);
      printf("Voff: %f\n",h->Voff);
      printf("bullock_spin: %f\n",h->bullock_spin);
      printf("halfmass_radius: %f\n",h->halfmass_radius);
#endif
   
#ifndef RUNSINGLE
// 44+5=49
if(ifile==0) printf("#time rhoh             x y z corevx corevy corevz  	    h_bulkvelx   h_bulkvely   h_bulkvelz 	    h_mass   h_r   h_rs   h_num_p   h_dens_tot   h_Xoff   h_Voff   h_kin_to_pot 	    h_energy   h_halfmass_radius 	    h_vrms   h_vmax   h_rvmax   h_spin   h_bullock_spin 	    h_Jx   h_Jy   h_Jz h_jx   h_jy   h_jz 	    h_jinx   h_jiny   h_jinz     h_djOverdTx   h_djOverdTy   h_djOverdTz 	    h_deOverdT 	    h_Ax   h_Ay   h_Az 	    h_b_to_a   h_c_to_a             h_sigma0 h_csigma h_sigmar h_sigmat h_sigmap\n");
      printf("%f %f %f %f %f %f %f %f %f %f %f %e %f %f %"PRId64" %"PRId64" %f %f %f %e %f %f %f %f %f %f %e %e %e %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
            0.1*((float)ifile),rhoh,         
	    cen[0],cen[1],cen[2],cen[3],cen[4],cen[5],         
	    h->bulkvel[0],h->bulkvel[1],h->bulkvel[2],         // 11
	    h->m,h->r,h->rs,h->num_p,h->dens_tot,h->Xoff,h->Voff,h->kin_to_pot,         // 19
	    h->energy,h->halfmass_radius,         //21
	    h->vrms,h->vmax,h->rvmax,h->spin,h->bullock_spin,         //26
	    h->J[0],h->J[1],h->J[2],         //29
	    h->J[0]/h->m,h->J[1]/h->m,h->J[2]/h->m,//32         
	    h->jin[0],h->jin[1],h->jin[2],         //35
	    h->djOverdT[0],h->djOverdT[1],h->djOverdT[2],        //38
	    h->deOverdT,         //39
	    h->A[0],h->A[1],h->A[2],         h->b_to_a,h->c_to_a,     //44
	    h->sigma0,h->csigma,sigmaH[0],sigmaH[1],sigmaH[2]);//49
#endif

//#################################################################################################################################
//#################################################################################################################################

#ifdef RUNSTARS 
      if(num_ps > 0){ // 2024.08.14  ZXY: if no star involved, do not run it to avoid crush ...... (2/3)

#ifdef DYDEBUG
      printf("#####################################################\n");
      printf("------------------ STAR ANALYSIS --------------------\n");
#endif
      // set number, particles, potentials
      total_p=num_ps;

      pot=pots;
      star = check_realloc(star,sizeof(struct halo),"Star");
      memset(star,0,sizeof(struct halo));
      memcpy(star->pos, cen, sizeof(float)*6);
      star->corevel[0]=cen[3]; star->corevel[1]=cen[4]; star->corevel[2]=cen[5];
      star->num_p=num_ps;
      //star->r=h->r;
   
      // Remove unbinding particles
      //qsort(pot, total_p, sizeof(struct potential), dist_compare_star);
      
      // 2024.05.28 ZXY: Now also remove unbound stars ......
      cenb[0]=0; cenb[1]=0; cenb[2]=0;
      cenunb[0]=0; cenunb[1]=0; cenunb[2]=0;
      nb=0; nunb=0;
      for (j=0; j<total_p; j++) {
        if ( (pot[j].ke-pot[j].pe) > 0 ) { // unbound particles
          nunb++;
          cenunb[0]+=pot[j].pos[0];
          cenunb[1]+=pot[j].pos[1];
          cenunb[2]+=pot[j].pos[2];
          total_p--;
          pot[j] = pot[total_p];
          j--;
        }
        else{ // bounded particles
           nb++;
           cenb[0]+=pot[j].pos[0];
           cenb[1]+=pot[j].pos[1];
           cenb[2]+=pot[j].pos[2];
           r2 = 0;
           for(k=0; k<3; k++) { dx = pot[j].pos[k]-cen[k]; r2 += dx*dx; }
        }
      }

//      nuseful=0;

      float Vprofsr[nticks-1];
      float Vprofst[nticks-1];
      float Vprofsp[nticks-1];
      float Sprofsr[nticks-1];
      float Sprofst[nticks-1];
      float Sprofsp[nticks-1];


//    ZXY 2022.07.03: Consider only bound star particles ......
      float rmax2, rmin2;
      rmax2 = 2;
      rmin2 = 0.005;

      nticks = 10;

      rx = rmin2;
      for(int ni=1;ni<nticks;ni++){ rx=rx*pow(rmax2/rmin2,1./((float)(nticks-1))); rpf2[ni]=rx; }

      getprof(&pot,total_p ,cenp,profs,rpf2,nticks);
      getMprof(&pot,total_p ,cenp,Mprofs,rpf2,nticks);
      getsigmaprof(&pot, total_p, cenp,Sprofsr,Sprofst,Sprofsp,rpf2,nticks);
      getVprof(&pot,total_p ,cenp,Vprofsr,Vprofst,Vprofsp,rpf2,nticks);

      char outpfs[1000];
      sprintf(outpfs,"%s/dataProfile/prof_star_%s_snap%.3d.txt", argv[1] ,argv[2],ifile);

      FILE *fstar = fopen(outpfs,"w");
      fprintf(fstar,"# r    rho    mr   sigmar   sigmat   sigmap   vr   vt   vp   beta\n");
      for(int i=0;i<nticks-1;i++) {
      fprintf(fstar,"%f    %f    %f    %f    %f    %f    %f    %f    %f\n",(rpf2[i]+rpf2[i+1])/2.,profs[i],Mprofs[i],Sprofsr[i],Sprofst[i],Sprofsp[i],Vprofsr[i],Vprofst[i],Vprofsp[i]);
      }
      fclose(fstar);

      // ZXY 2022/01/25: In pots, r2 has already been updated using the potential minimum as center ......
      qsort(pot, total_p, sizeof(struct potential), dist_compare_sg);
      star->num_p=total_p; // update star particle number

#ifdef DYDEBUG
   
      printf("-------------------------------------------------------\n");
      printf("Number total_p %d\n",total_p);
      //printf("Number of nb: %d total_p %d\n",nb,total_p);
      //printf("Binding particle fraction:           %f\n",(float)total_p/((float)num_ps));
      //printf("Fraction of star particles used:     %f\n",(float)total_p/((float)num_ps));
      //printf(">>>> Mass center of binding particles:   %f %f %f\n",cenb[0]/((float)total_p),cenb[1]/((float)total_p),cenb[2]/((float)total_p));
      //printf(">>>> Mass center of unbinding particles: %f %f %f\n",cenunb[0]/((float)(num_ps-total_p)),cenunb[1]/((float)(num_ps-total_p)),cenunb[2]/((float)(num_ps-total_p)));
#endif   
      calc_mass_definition_star();
      calc_basic_halo_props_star(star);
      calc_additional_halo_props_star(star);

      char stec[2000];
      sprintf(stec,"%s/dataProfile/star_eccentricity_%.3d.txt", argv[1], isnap);
      printf("\n\nStar eccentricity profile recorded in %s \n\n",stec);
      FILE *fse = fopen(stec,"w");
      for(i=0; star->starec[0][i] != 0; i++) {
        fprintf(fse,"%f  %f  %f  %f\n", 1000*star->starec[0][i], star->starec[1][i], star->starec[2][i], star->starec[3][i]);
      }
      fclose(fse);

      // ZXY 2022.01.29: After calc_additional_halo_props, we get bound mass (halo + star) within rvir 
      //                 (h->mgrav and star->mgrav) ......
      // --------------------------------------------------------------------------------------------
      //
/*
      char bdmass[1000];

      FILE *boundmass;

      sprintf(bdmass,"%s/boundmass.txt",argv[1]);
      boundmass = fopen(bdmass,"a+");

      r2 = sqrt( cen[0]*cen[0] + cen[1]*cen[1] + cen[2]*cen[2] );
      float r2dm = sqrt(cendm[0]*cendm[0]+cendm[1]*cendm[1]+cendm[2]*cendm[2]);

      fprintf(boundmass, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", time_snap, h->mgrav, mdm_Rh, cen[0], cen[1], cen[2], r2, SELF_PE, SELF_KE, TOTAL_PE, TOTAL_KE, mdm_10, mdm_20, mdm_50, BOUND_PE, BOUND_KE);

      fclose(boundmass);
      printf("ZXY: bound mass recorded ..........\n\n");
      // --------------------------------------------------------------------------------------------
*/
      // Note: to avoid divergence, if alpha=-1, r is softed...
      float mtme0,mtmem1,mtmep1;
      float sigmarp;
      sigmarp=getv0(&pot,total_p,cen,0.0,7.6,-1);
      mtmem1=0.572958*1e10*(sigmarp*sigmarp)*7.6*7.6/43009.3;
      sigmarp=getv0(&pot,total_p,cen,0.0,7.6,0);
      mtme0 =1.9     *1e10*(sigmarp*sigmarp)*7.6/43009.3;
      sigmarp=getv0(&pot,total_p,cen,0.0,7.6,1);
      mtmep1=3.69239*1e10*(sigmarp*sigmarp)*1/43009.3;

      float mptme1;
      sigmarp=getv0(&pot,total_p,cen,0.0,7.0,1);
      mptme1 =6.3662     *1e10*(sigmarp*sigmarp)*1/43009.3;

      float mh70=PARTICLE_MASS*(nh70);
      float ms70=PARTICLE_MASS*(ns70);
      float mtot70=PARTICLE_MASS*(nh70+ns70);
      float mh76=PARTICLE_MASS*(nh76);
      float ms76=PARTICLE_MASS*(ns76);
      float mtot76=PARTICLE_MASS*(nh76+ns76);
      
      float sigmaS[3]={0,0,0};
      getsigma(&pot,total_p,cen,sigmaS,0.0, 7 /*  star->r  */);


      // ZXY 2022.06.29: Calculate stellar velocity dispersion in anothor way ......
/*
      float vsigma[nticks-1];
      for(int i=0; i<nticks-1; i++){
         vsigma[i] = getsigma(&pot,total_p,cen,sigmaH,rpf2[i],rpf2[i+1]);
         printf("ZXY: vsigma in %f to %f kpc: %f ...\n", rpf2[i], rpf2[i+1], vsigma[i]);
      }
      char outpvs2[1000];
      sprintf(outpvs2,"%s/dataProfile/prof_star_vdisp_%s_%.3d.txt", argv[1], argv[2], ifile);
      FILE *fvs   = fopen(outpvs2,"w");

      fprintf(fvs,"# r    rho   mr   vdisp\n");
      for(int i=0;i<nticks-1;i++) {
        fprintf(fvs,  "%f    %f    %f   %f\n",(rpf[i]+rpf[i+1])/2.,profh[i],Mprofh[i],vsigma[i]);
      }
      fclose(fvs);
*/


#ifdef DYDEBUG
      printf("-------------------------------------------------------\n");
      printf("\n----RockStar results for the stellar component ------\n");
      printf("pos[6]: %f %f %f %f %f %f\n",star->pos[0],star->pos[1],star->pos[2],star->pos[3],star->pos[4],star->pos[5]);
      printf("core vel %f %f %f \n",star->corevel[0],star->corevel[1],star->corevel[2]);
      printf("bulk vel %f %f %f \n",star->bulkvel[0],star->bulkvel[1],star->bulkvel[2]);
      printf("mass: %e\n",star->m); // vir
      printf("mbound: %e\n",star->mgrav);
      printf("r: %f\n",star->r); // vir
      printf("center: %f, %f, %f \n", star->pos[0], star->pos[1], star->pos[2]);
      printf("rs: %f\n",star->rs);
      printf("num_p: %d\n",star->num_p);
      printf("dens_tot: %d\n",star->dens_tot);
      printf("vrms: %f\n",star->vrms);
      printf("vmax: %f\n",star->vmax);
      printf("rvmax: %f\n",star->rvmax);
      printf("J[3] %e %e %e\n",star->J[0],star->J[1],star->J[2]);
      printf("j[3] %f %f %f\n",star->J[0]/star->m,star->J[1]/star->m,star->J[2]/star->m);
      printf("jin[3] (r<7.6 kpc) %f %f %f\n",star->jin[0],star->jin[1],star->jin[2]);
      printf("djOverdT[3] %f %f %f\n",star->djOverdT[0],star->djOverdT[1],star->djOverdT[2]);
      printf("deOverdT: %f\n",star->deOverdT);
      printf("b_to_a: %f\n",star->b_to_a);
      printf("c_to_a: %f\n",star->c_to_a);
      printf("A[3]: %f %f %f\n",star->A[0],star->A[1],star->A[2]);
      printf("kin_to_pot: %f\n",star->kin_to_pot);
      printf("energy: %e\n",star->energy);
      printf("spin: %f\n",star->spin);
      printf("Xoff: %f\n",star->Xoff);
      printf("Voff: %f\n",star->Voff);
      printf("bullock_spin: %f\n",star->bullock_spin);
      printf("halfmass_radius: %f\n",star->halfmass_radius);
#endif
      float vhx=h->A[0];    float vhy=h->A[1];
      float vsx=star->A[0]; float vsy=star->A[1];
      float phase = acos((vhx*vsx+vhy*vsy)/(sqrt(vhx*vhx+vhy*vhy)*sqrt(vsx*vsx+vsy*vsy)));


      //if(ifile==0) printf("#time rhoh x y hAx hAy sAx sAy phase\n");
      //printf("%f %f %f %f %f %f %f %f %f\n",0.1*((float)ifile),rhoh,cen[0],cen[1],h->A[0],h->A[1],star->A[0],star->A[1],phase);
      // Almost full information...
/*
      if(ifile==0) printf("#time rhoh phase
           x y z corevx corevy corevz
           h_bulkvelx   h_bulkvely   h_bulkvelz 
           h_mass   h_r   h_rs   h_num_p   h_dens_tot   h_Xoff   h_Voff   h_kin_to_pot 
           h_energy   h_halfmass_radius 
           h_vrms   h_vmax   h_rvmax   h_spin   h_bullock_spin 
           h_Jx   h_Jy   h_Jz 
           h_jx   h_jy   h_jz 
           h_jinx   h_jiny   h_jinz 
           h_djOverdTx   h_djOverdTy   h_djOverdTz 
           h_deOverdT 
           h_Ax   h_Ay   h_Az 
           h_b_to_a   h_c_to_a 
           s_bulkvelx   s_bulkvely   s_bulkvelz 
           s_mass   s_r   s_rs   s_num_p   s_dens_tot   s_Xoff   s_Voff   s_kin_to_pot 
           s_energy   s_halfmass_radius 
           s_vrms   s_vmax   s_rvmax   s_spin   s_bullock_spin 
           s_Jx   s_Jy   s_Jz 
           s_jx   s_jy   s_jz 
           s_jinx   s_jiny   s_jinz 
           s_djOverdTx   s_djOverdTy   s_djOverdTz 
           s_deOverdT 
           s_Ax   s_Ay   s_Az 
           s_b_to_a   s_c_to_a \n");
        printf("%f %f %f, %f %f %f %f %f %f, 
%f %f %f, %e %f %f %"PRId64" %"PRId64" %f %f %f, %e %f, %f %f %f %f %f, %e %e %e, %f %f %f, %f %f %f, %f %f %f, %f, %f %f %f, %f %f, 
%f %f %f, %e %f %f %"PRId64" %"PRId64" %f %f %f, %e %f, %f %f %f %f %f, %e %e %e, %f %f %f, %f %f %f, %f %f %f, %f, %f %f %f, %f %f, 
\n",*/
      float dist0=sqrt(cen[0]*cen[0]+cen[1]*cen[1]+cen[2]*cen[2])*1e3;
      float vel0=sqrt(cen[3]*cen[3]+cen[4]*cen[4]+cen[5]*cen[5]);
#ifndef RUNSINGLE
if(ifile==0) printf("#time rhoh rhos phase x y z corevx corevy corevz h_bulkvelx   h_bulkvely   h_bulkvelz h_mass   h_r   h_rs   h_num_p   h_dens_tot   h_Xoff   h_Voff   h_kin_to_pot h_energy   h_halfmass_radius h_vrms   h_vmax   h_rvmax   h_spin   h_bullock_spin h_Jx   h_Jy   h_Jz h_jx   h_jy   h_jz h_jinx   h_jiny   h_jinz h_djOverdTx   h_djOverdTy   h_djOverdTz h_deOverdT h_Ax   h_Ay   h_Az h_b_to_a   h_c_to_a s_bulkvelx   s_bulkvely   s_bulkvelz s_mass   s_r   s_rs   s_num_p   s_dens_tot   s_Xoff   s_Voff   s_kin_to_pot s_energy   s_halfmass_radius s_vrms   s_vmax   s_rvmax   s_spin   s_bullock_spin s_Jx   s_Jy   s_Jz s_jx   s_jy   s_jz s_jinx   s_jiny   s_jinz s_djOverdTx   s_djOverdTy   s_djOverdTz s_deOverdT s_Ax   s_Ay   s_Az s_b_to_a   s_c_to_a dist vel mtmem1 mtme0 mtmep1 mh76 ms76 mtot76 h_sigma0 h_csigma h_sigmar h_sigmat h_sigmap s_sigmar s_sigmat s_sigmap\n");
      printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %e %f %f %"PRId64" %"PRId64" %f %f %f %e %f %f %f %f %f %f %e %e %e %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %e %f %f %"PRId64" %"PRId64" %f %f %f %e %f %f %f %f %f %f %e %e %e %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",snapshot2time*((float)ifile),rhoh,rhos,phase,         cen[0],cen[1],cen[2],cen[3],cen[4],cen[5],         h->bulkvel[0],h->bulkvel[1],h->bulkvel[2],         h->m,h->r,h->rs,h->num_p,h->dens_tot,h->Xoff,h->Voff,h->kin_to_pot,         h->energy,h->halfmass_radius,         h->vrms,h->vmax,h->rvmax,h->spin,h->bullock_spin,         h->J[0],h->J[1],h->J[2],         h->J[0]/h->m,h->J[1]/h->m,h->J[2]/h->m,         h->jin[0],h->jin[1],h->jin[2],         h->djOverdT[0],h->djOverdT[1],h->djOverdT[2],         h->deOverdT,         h->A[0],h->A[1],h->A[2],         h->b_to_a,h->c_to_a,         star->bulkvel[0],star->bulkvel[1],star->bulkvel[2],         star->m,star->r,star->rs,star->num_p,star->dens_tot,star->Xoff,star->Voff,star->kin_to_pot,         star->energy,star->halfmass_radius,         star->vrms,star->vmax,star->rvmax,star->spin,star->bullock_spin,         star->J[0],star->J[1],star->J[2],         star->J[0]/star->m,star->J[1]/star->m,star->J[2]/star->m,         star->jin[0],star->jin[1],star->jin[2],         star->djOverdT[0],star->djOverdT[1],star->djOverdT[2],         star->deOverdT,         star->A[0],star->A[1],star->A[2],         star->b_to_a,star->c_to_a,dist0,vel0,mtmem1,mtme0,mtmep1,mh76,ms76,mtot76,h->sigma0,h->csigma,sigmaH[0],sigmaH[1],sigmaH[2],sigmaS[0],sigmaS[1],sigmaS[2]);
#endif

#ifdef RUNSINGLE
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("XXXXXXX halo mass %e, mh70= %e, mh76= %e, rho0=%f, rH = %f, r = %f, rs = %f, b_to_a = %f, c_to_a = %f, A = %f\n",h->m,mh70,mh76,rhoh,h->halfmass_radius,h->r,h->rs,h->b_to_a,h->c_to_a,sqrt(h->A[0]*h->A[0]+h->A[1]*h->A[1]+h->A[2]*h->A[2]));
      printf("XXXXXXX star mass %e, ms76= %e, rho0=%f, rH = %f, r = %f, rs = %f, b_to_a = %f, c_to_a = %f, A = %f\n",star->m,ms76,rhos,star->halfmass_radius,star->r,star->rs,star->b_to_a,star->c_to_a,sqrt(star->A[0]*star->A[0]+star->A[1]*star->A[1]+star->A[2]*star->A[2]));
//      printf("XXXXXXX check sigmarp: %f, TME mass: %e, mtot7.6: %e\n",sigmarp,star->mtme,mtot76);
      printf("XXXXXXX t = %f, rH = %f, ms76 = %e, mtmem1 = %e, mtme0 = %e, mtmep1 = %e, mtot76 = %e, (80/d)^2+(293/v)^2 = %f \n",0.1*((float)ifile)*0.5,star->halfmass_radius,ms76,mtmem1,mtme0,mtmep1,mtot76,(80./dist0)*(80./dist0)+(293./vel0)*(293./vel0));
      //printf(" >>>>>>> rH_star = %f, b/a = %f \n",star->halfmass_radius,star->b_to_a);
      //printf("DF-2 like: \n");
      printf("    %1.1f      &       %1.2f          &   %1.2f      &             %1.2f                  &        %1.2f                              &         %1.2f                          &    %1.1f      \n",0.1*((float)ifile)*0.5,star->halfmass_radius,star->b_to_a,ms76/1e8,mtot76/1e8,mtme0/1e8,(80./dist0)*(80./dist0)+(293./vel0)*(293./vel0));
      //printf("New TME mass: %e\n",mtme0);
      printf("New bound halo mass: %e, stellar mass: %e, TME mass: %e, r: %f\n",h->m,star->m,mtme0,h->r);
      printf("XXXXXXX halo x=%f, y=%f, z=%f, vx=%f, vy=%f, vz=%f\n",cen[0]*1000,cen[1]*1000,cen[2]*1000,cen[3],cen[4],cen[5]);
      printf("PARTICLE_MASS %f ns76 %f\n",PARTICLE_MASS,ns76);
      //printf("XXXXXXX check Mass, sigma0 = %f, Mtot = %e, Vcirc = %f\n",h->sigma0,3.5e11*pow(h->sigma0/220,4.43),sqrt(2.)*h->sigma0);
      //printf("XXXXXXX check t = %f, M200h = %e Msun, Vmax = %f, rvmax = %f, 2.16*rs = %f, sigma0=%f \n",0.1*((float)ifile)*0.5,h->m,h->vmax,h->rvmax,2.16*h->rs,h->sigma0);
      printf("M200 = %e Msun, r200 = %f, rs = %f, \n\t Vmax = %f, rvmax = %f, 2.16*rs = %f, rx=%f, vy=%f \n",h->m,h->r,h->rs,h->vmax,h->rvmax,2.16*h->rs,cen[0],cen[4]);
      printf(" vr=%f, vth=%f, vphi=%f \n ",sigmaH[0],sigmaH[1],sigmaH[2]);
      //printf(" CHECK if ms76*2 %e > mtme0 %e, vel = %f\n",ms76*2,mtme0,vel0);
      //printf("DF-4 like: \n");
      //printf("    %1.1f      &             %1.2f                  &        %1.2f                              &         %1.2f                          &    %1.1f    &  O  \n",0.1*((float)ifile)*0.5,ms70/1e8,mtot70/1e8,mptme1/1e8,vel0);
      //printf(" CHECK if ms70*2 %e > mptme1 %e, vel = %f\n",ms70*2,mptme1,vel0);
      //printf("        cx = %f; cy = %f; cz = %f; cvx = %f; cvy=%f; cvz = %f;\n",cen[0]*1000,cen[1]*1000,cen[2]*1000,cen[3],cen[4],cen[5]);
      //printf("        dist = %f, velocity = %f, (80/d)^2+(293/v)^2 = %f (?= 1)\n",dist0,vel0,(80./dist0)*(80./dist0)+(293./vel0)*(293./vel0));
#endif

  } // End of num_ps > 0 ...... (2/3)

#endif // end RUNSTARS

      char bdmass[1000];

      FILE *boundmass;

      sprintf(bdmass,"%s/boundmass.txt",argv[1]);
      boundmass = fopen(bdmass,"a+");

      float r2sq = sqrt( cen[0]*cen[0] + cen[1]*cen[1] + cen[2]*cen[2] );

      // todo
       /* time_snap: (the snapshot_#) * 0.2
        * h->mgrav: bound mass ?
        * mdm_Rh:
        * cen[0-2]: position with potential minimum; cen[3-5]: velocity with potential minimum ?
        * r2sq: distance^2
        * SELF_PE: potential ?
        * SELF_KE: kinetic energy ?
        * TOTAL_PE: SELF_PE + EXT_PE
        * TOTAL_KE: sum over kinetic energy of particles without far away with the halo.
        *           TOTAL_KE = relative kinetic energy based on center point + absolute kinetic energy ?
        * mdm_x: the total masses of particles whose distance < (x pc)
        * BOUND_PE: bound potential = self potential ?
        * BOUND_KE: bound kinetic energy = self kinetic energy ?
        * */
      fprintf(boundmass, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
              time_snap, h->mgrav, mdm_Rh, cen[0], cen[1], cen[2], r2sq, SELF_PE, SELF_KE, TOTAL_PE, TOTAL_KE, mdm_10,
              mdm_20, mdm_50, BOUND_PE, BOUND_KE);

      fclose(boundmass);
      printf("ZXY: bound mass recorded ..........\n\n");

//      Hz=h0/9.777752*sqrt(1.0-Om + Om/pow(SCALE_NOW,3));
//      dtime=(SCALE_NOW-scale_m1)/(Hz*pow(SCALE_NOW,3)); 
//      scale_m1 = SCALE_NOW;

      //cen[0]-=cen[3]*SCALE_NOW*dtime*1e-3;
      //cen[1]-=cen[4]*SCALE_NOW*dtime*1e-3;
      //cen[2]-=cen[5]*SCALE_NOW*dtime*1e-3;

      free(pdm);
      printf("t1 ..........................\n");
      free(ps);
      printf("t2 ..........................\n");
      free(h);
      printf("t3 ..........................\n");
//      free(poth);
//      printf("t4 ..........................\n");
#ifdef RUNSTARS
  if(num_ps > 0)  // 2024.08.14  ZXY: if no star involved, do not run it to avoid crush ...... (3/3)
      free(star);
#endif
      free(pots);
      printf("t5 ..........................\n");
//     free(pot);
//     printf("t6 ..........................\n");
      free(p);
//      p=NULL;
//      pot=NULL;
      printf("End ..........................\n");
   }
   
   return 0;
}
