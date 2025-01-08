#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "io/io_gadget.h"
#include "io/getcenter.h"
#include "io/io_generic.h"
#include "particle.h"
#include "check_syscalls.h"
#include "myutils.h"
#include "potential.h"
#include "groupies.h"
#include "universal_constants.h"

#include "singleprop.c"

int runsingle(, char *filename){

   load_particles_gadget2(filename, &pdm, &num_ph);
   load_particles_gadgetStar2(filename, &ps, &num_ps);
   printf("N halos: %d, N stars: %d\n",num_ph,num_ps);
   //for(i=0;i<1;i++){
   //   printf("Check id pos for the 1st particle:\n \t %d %f %f %f %f %f %f\n",p[i].id,p[i].pos[0],p[i].pos[1],p[i].pos[2],p[i].pos[3],p[i].pos[4],p[i].pos[5]);
   //}


   p = (struct particle *) check_realloc(p, (num_ph+num_ps)*sizeof(struct particle), "Allocating room for particles.");

   for(i=0;i<num_ph;i++){ p[i].id=i; memcpy(p[i].pos, pdm[i].pos, sizeof(float)*6); }
   for(j=0;j<num_ps;j++){ p[num_ph+j].id=num_ph+j; memcpy(p[num_ph+j].pos, ps[j].pos, sizeof(float)*6); }
   //printf("XXXXX %f %f %f %f",p[0].pos[0],p[num_ph-1].pos[0],p[num_ph].pos[0],p[num_ph+num_ps-1].pos[0]);
   
   int64_t total_p = num_ph+num_ps;
   float cen[6];
   // calc center position and velocity
   //getcenter(&pdm,num_ph,cen,0.6);
   //getcenter(&ps,num_ps,cen,0.6);
   double rho0 = getcenter(&p,total_p,cen,0.6);
   printf(">>>>>> Central core density: %g\n",rho0);

   // We calculate potential only once, using the center just found. 
   float dx, r2;
   pot = (struct potential *) check_realloc(pot, (total_p)*sizeof(struct potential), "Allocating room for potentials.");
   poth = (struct potential *) check_realloc(poth, (num_ph)*sizeof(struct potential), "Allocating room for potentials.");
   pots = (struct potential *) check_realloc(pots, (num_ps)*sizeof(struct potential), "Allocating room for potentials.");
   memset(pot, 0, sizeof(struct potential)*(total_p));
   memset(poth, 0, sizeof(struct potential)*(num_ph));
   memset(pots, 0, sizeof(struct potential)*(num_ps));

   for(j=0; j<total_p; j++) {
     if(j<num_ph) pot[j].type=1;
     else pot[j].type=4;
     r2 = 0;
     for (k=0; k<3; k++) { dx=p[j].pos[k] - cen[k]; r2+=dx*dx; }
     pot[j].r2 = r2;
     if(r2==0) printf("Warning!!!!! r2==0");
     memcpy(pot[j].pos, p[j].pos, sizeof(float)*6);
   }
   qsort(pot, total_p, sizeof(struct potential), dist_compare_sg);

   //Assumes center + velocity already calculated.
   setMOR(cen);
   compute_potential(pot,total_p);
   float corevel[3]={cen[3],cen[4],cen[5]};
   compute_kinetic_energy(pot, total_p, corevel, cen);

   j=0;k=0;
   for(i=0; i<total_p; i++) {
      if(pot[i].type==1){ memcpy(&poth[j],&pot[i],sizeof(struct potential)*1); j++;}
      if(pot[i].type==4){ memcpy(&pots[k],&pot[i],sizeof(struct potential)*1); k++;}
   }
   //printf("num_ph == %d; num_ps == %d \n",j,k);
   free(pot);

   // SIDM 101: consider only the inner 20 kpc...
   //--------------------------------------------------------

   printf("------------------ HALO ANALYSIS --------------------\n");

   //---------------------------------------------
   //struct halo *h = malloc(sizeof(struct halo));
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

   FILE *fb   = fopen("outbind.txt","w");
   FILE *funb = fopen("outunbind.txt","w");
   //qsort(pot, total_p, sizeof(struct potential), dist_compare_sg);
   //calculate_corevel_sg(h, pot, total_p);

#ifdef STRIPPED
   float rbound=20e-3;
   printf(">>>>>>>> We consider only particles in a sphere of r = %f around the halo\n",rbound);
#endif
   for (j=0; j<total_p; j++) {
#ifdef STRIPPED
     if ( (pot[j].ke-pot[j].pe) > 0 || pot[j].r2>rbound*rbound) {
#else
     if ( (pot[j].ke-pot[j].pe) > 0 ) {
#endif
       if((pot[j].ke-pot[j].pe) > 0){ 
          nunb++; 
          if(nunb%1000==0) fprintf(funb,"%f %f %f\n",pot[j].pos[0],pot[j].pos[1],1.0);
          cenunb[0]+=pot[j].pos[0];
          cenunb[1]+=pot[j].pos[1];
          cenunb[2]+=pot[j].pos[2];
       }
       if((pot[j].ke-pot[j].pe) <= 0){ 
          nb++;
          if(nb%1000==0) fprintf(fb,"%f %f %f\n",pot[j].pos[0],pot[j].pos[1],1.0);
          cenb[0]+=pot[j].pos[0];
          cenb[1]+=pot[j].pos[1];
          cenb[2]+=pot[j].pos[2];
       }
       total_p--;
       pot[j] = pot[total_p];
       j--;
     }
     else{
        nb++;
        if(nb%1000==0) fprintf(fb,"%f %f %f\n",pot[j].pos[0],pot[j].pos[1],1.0);
        cenb[0]+=pot[j].pos[0];
        cenb[1]+=pot[j].pos[1];
        cenb[2]+=pot[j].pos[2];
     }
   }
   fclose(fb);
   fclose(funb);
   qsort(pot, total_p, sizeof(struct potential), dist_compare_sg);
   h->num_p=total_p; // update halo particle number

   printf("-------------------------------------------------------\n");
   printf("Number of nb: %d\n",nb);
   printf("Binding particle fraction:               %f\n",(float)nb/((float)num_ph));
   printf("Halo particle fraction:                  %f\n",(float)total_p/((float)num_ph));
   printf(">>>> Mass center of binding particles:   %f %f %f\n",cenb[0]/((float)total_p),cenb[1]/((float)total_p),cenb[2]/((float)total_p));
   printf(">>>> Mass center of unbinding particles: %f %f %f\n",cenunb[0]/((float)(num_ph-total_p)),cenunb[1]/((float)(num_ph-total_p)),cenunb[2]/((float)(num_ph-total_p)));

   //printf("Halo center: x=%f, y=%f, z=%f \n",h[0].pos[0],h->pos[1],h->pos[2]);
   //for(j=0; j<20; j++) {
   //   printf("check %f %f\n",pot[j].pos[0],pot[j].pos[4]);
   //}
   //return 0;

   calc_mass_definition_sg();
   calc_basic_halo_props_sg(h);
   calc_additional_halo_props_sg(h);

   //printf("TEST main: %f\n",arrMOverR[100]);
   //compute_potential(pot,10000);

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
   printf("c_to_a: %f\n",h->b_to_a);
   printf("A[3]: %f %f %f\n",h->A[0],h->A[1],h->A[2]);
   printf("kin_to_pot: %f\n",h->kin_to_pot);
   printf("energy: %e\n",h->energy);
   printf("spin: %f\n",h->spin);
   printf("Xoff: %f\n",h->Xoff);
   printf("Voff: %f\n",h->Voff);
   printf("bullock_spin: %f\n",h->bullock_spin);
   printf("halfmass_radius: %f\n",h->halfmass_radius);

   printf("#####################################################\n");
   printf("------------------ STAR ANALYSIS --------------------\n");
   // set number, particles, potentials
   total_p=num_ps;
   pot=pots;

   star = check_realloc(star,sizeof(struct halo),"Star");
   memset(star,0,sizeof(struct halo));
   memcpy(star->pos, cen, sizeof(float)*6);
   star->corevel[0]=cen[3]; star->corevel[1]=cen[4]; star->corevel[2]=cen[5];
   star->num_p=num_ps;
   star->r=h->r;

   // Remove unbinding particles
   //qsort(pot, total_p, sizeof(struct potential), dist_compare_sg);
   cenb[0]=0;cenb[1]=0;cenb[2]=0;
   cenunb[0]=0;cenunb[1]=0;cenunb[2]=0;
   nb=0,nunb=0;
   for (j=0; j<total_p; j++) {
     if ( (pot[j].ke-pot[j].pe) > 0) {
       cenunb[0]+=pot[j].pos[0];
       cenunb[1]+=pot[j].pos[1];
       cenunb[2]+=pot[j].pos[2];
       total_p--;nunb++;
       pot[j] = pot[total_p];
       j--;
     }
     else{
        nb++;
        cenb[0]+=pot[j].pos[0];
        cenb[1]+=pot[j].pos[1];
        cenb[2]+=pot[j].pos[2];
     }
   }
   qsort(pot, total_p, sizeof(struct potential), dist_compare_sg);
   star->num_p=total_p; // update star particle number

   printf("-------------------------------------------------------\n");
   printf("Number of nb: %d total_p %d\n",nb,total_p);
   printf("Binding particle fraction:           %f\n",(float)total_p/((float)num_ps));
   printf("Fraction of star particles used:     %f\n",(float)total_p/((float)num_ps));
   printf(">>>> Mass center of binding particles:   %f %f %f\n",cenb[0]/((float)total_p),cenb[1]/((float)total_p),cenb[2]/((float)total_p));
   printf(">>>> Mass center of unbinding particles: %f %f %f\n",cenunb[0]/((float)(num_ps-total_p)),cenunb[1]/((float)(num_ps-total_p)),cenunb[2]/((float)(num_ps-total_p)));

   //calc_mass_definition_sg();
   calc_basic_halo_props_sg(star);
   calc_additional_halo_props_sg(star);


   printf("-------------------------------------------------------\n");
   printf("\n----RockStar results for the stellar component ------\n");
   printf("pos[6]: %f %f %f %f %f %f\n",star->pos[0],star->pos[1],star->pos[2],star->pos[3],star->pos[4],star->pos[5]);
   printf("core vel %f %f %f \n",star->corevel[0],star->corevel[1],star->corevel[2]);
   printf("bulk vel %f %f %f \n",star->bulkvel[0],star->bulkvel[1],star->bulkvel[2]);
   printf("mass: %e\n",star->m); // vir
   printf("mbound: %e\n",star->mgrav);
   printf("r: %f\n",star->r); // vir
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
   printf("c_to_a: %f\n",star->b_to_a);
   printf("A[3]: %f %f %f\n",star->A[0],star->A[1],star->A[2]);
   printf("kin_to_pot: %f\n",star->kin_to_pot);
   printf("energy: %e\n",star->energy);
   printf("spin: %f\n",star->spin);
   printf("Xoff: %f\n",star->Xoff);
   printf("Voff: %f\n",star->Voff);
   printf("bullock_spin: %f\n",star->bullock_spin);
   printf("halfmass_radius: %f\n",star->halfmass_radius);

   free(p);
   free(pdm);
   free(ps);
   free(h);
   free(star);
   free(poth);
   free(pots);
   //free(pot);
   return 0;
}
