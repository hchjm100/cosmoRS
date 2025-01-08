#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../universal_constants.h"
#include "../config_vars.h"
#include "getrho.h"

float getrho(struct particle **p,int64_t num_p,float *cen,float rmax){

    // Convert LU into kpc, results in Gadget units: 1e10 Msun/kpc^3

    double x,y,z;
    double vx,vy,vz;
    double vx1,vy1,vz1;
    const double rin=0.1;
    const double frac=0.95;
    float mass = PARTICLE_MASS/1e10;

    double npart=0;
    double cx0,cy0,cz0;
    double cx1,cy1,cz1;
    double rden,rho0;
    int64_t i;

    const int nstep=20;
    double r=0;
    double sr[nstep];
    double sw[nstep];
    int n=0;

    cx0=cen[0]*1000; cy0=cen[1]*1000; cz0=cen[2]*1000; 
    sr[0]=0;
    r=rin;
    for(n=1;n<=nstep;n++){
       r=r*pow(rmax/rin,1./((double*)(nstep-1)));
       sr[n]=r;
       sw[n-1]=mass/(4./3.*M_PI*(pow(sr[n],3)-pow(sr[n-1],3)));
    }

    for(i=0; i<num_p ;i++){
        x=(*p)[i].pos[0];
        y=(*p)[i].pos[1];
        z=(*p)[i].pos[2];
        r2 = (x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0);
        r2 = r2*1000*1000*10*10;
        if(r2<1.0) npt[0]++;

    // =======================================================
    
    // Find a radius that contain ~1000 particles inside
    // Condition: 1 kpc or 1000 particles
    int64_t npt[20]={0};
    double r2=0;
    for(i=0; i<num_p ;i++){
        x=(*p)[i].pos[0];
        y=(*p)[i].pos[1];
        z=(*p)[i].pos[2];
        r2 = (x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0);
        r2 = r2*1000*1000*10*10;
        if(r2<1.0) npt[0]++;
        else if(r2<2*2) npt[1]++; else if(r2<3*3) npt[2]++; else if(r2<4*4) npt[3]++; else if(r2<5*5) npt[4]++; else if(r2<6*6) npt[5]++; 
        else if(r2<7*7) npt[6]++; else if(r2<8*8) npt[7]++; else if(r2<9*9) npt[8]++; else if(r2<10*10) npt[9]++;
        else if(r2<11*11) npt[10]++; else if(r2<12*12) npt[11]++; else if(r2<13*13) npt[12]++; else if(r2<14*14) npt[13]++; 
        else if(r2<15*15) npt[14]++; else if(r2<16*16) npt[15]++; else if(r2<17*17) npt[16]++; else if(r2<18*18) npt[17]++; 
        else if(r2<19*19) npt[18]++; else if(r2<20*20) npt[19]++;
    }
    double rmin=0.1; // kpc

    int64_t nmin=0;
    for(i=0;i<20;i++){
       if(nmin<200) nmin+=npt[i];
       else{ rmin=0.0001*(double)(i+1); break; }
       //printf("%d \n",npt[i]);
    }

    // A new iteration to obtain velocities
    npart=0;
    vx1=0;vy1=0;vz1=0;
    for(i=0; i<num_p ;i++){
	x=(*p)[i].pos[0];
	y=(*p)[i].pos[1];
	z=(*p)[i].pos[2];
	vx=(*p)[i].pos[3];
	vy=(*p)[i].pos[4];
	vz=(*p)[i].pos[5];
	if((x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0)<rmin*rmin){
	    vx1+=vx;
	    vy1+=vy;
	    vz1+=vz;
	    npart++;
	}
    }
    vx1=vx1/((double)(npart));
    vy1=vy1/((double)(npart));
    vz1=vz1/((double)(npart));
    rho0=PARTICLE_MASS*((double)npart)/(4./3. *M_PI*rmin*rmin*rmin*1e18);

    return rho0;
}
