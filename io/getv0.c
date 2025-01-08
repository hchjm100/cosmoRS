#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../universal_constants.h"
#include "../config_vars.h"
#include "getv0.h"

// 1D velocity dispersion at rin+rout
float getv0(struct potential **p,int64_t num_p,float *cen,float rin,float rout,float alpha){

    float x,y,z;
    float vx,vy,vz,vr,vt,vp;
    float vxbar=0,vybar=0,vzbar=0,vrbar=0,vtbar=0,vpbar=0;
    float sigma0=0,th=0,phi=0;

    float npart=0;
    float r2=0;
    float cx0,cy0,cz0;
    float vx0,vy0,vz0;
    int64_t i;
    int nxx=0;

    cx0=cen[0];   cy0=cen[1];   cz0=cen[2]; 
    vx0=cen[3];   vy0=cen[4];   vz0=cen[5];

    // particles sort by r2... 
    for(i=0; i<num_p ;i++){
        x=(*p)[i].pos[0];
        y=(*p)[i].pos[1];
        z=(*p)[i].pos[2];
        vx=(*p)[i].pos[3]-vx0;
        vy=(*p)[i].pos[4]-vy0;
        vz=(*p)[i].pos[5]-vz0;

        float r;
        if(alpha<0) r = sqrt((x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0)+0.05*0.05*1e-6); // Mpc
        else r = sqrt((x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0));
        if(r<0.001) nxx++;
        //float r = sqrt((*p)[i].r2+1e-10);
        //if((r*1000.-rin) > 0 && (r*1000.-(rin+rout))<0 && ((*p)[i].ke-(*p)[i].pe) < 0){
        if((r*1000.-rin) > 0 && (r*1000.-(rin+rout))<0){
           // error: need to calculate vx,vy,vz relative to the bulk motion... ... 
           // --------------------------------------------------------------------
           th=acos((z-cz0)/r);
           phi=atan2((y-cy0),(x-cx0));
           vr = vx*cos(phi)*sin(th)+vy*sin(phi)*sin(th)+vz*cos(th);
           vt = vx*cos(phi)*cos(th)+vy*sin(phi)*cos(th)-vz*sin(th);
           vp = -vx*sin(phi)+vy*cos(phi);

           vrbar+=vr;
           vtbar+=vt;
           vpbar+=vp;
           vxbar+=vx; 
           vybar+=vy; 
           vzbar+=vz; 
           npart++;
        }
    }
    //printf("xxxxxxxxxxxx inner 1kpc has %d particles\n",nxx);
    vxbar=vxbar/npart;
    vybar=vybar/npart;
    vzbar=vzbar/npart;
    vrbar=vrbar/npart;
    vtbar=vtbar/npart;
    vpbar=vpbar/npart;

    for(i=0; i<num_p ;i++){
	x=(*p)[i].pos[0];
	y=(*p)[i].pos[1];
	z=(*p)[i].pos[2];
        vx=(*p)[i].pos[3]-vx0;
        vy=(*p)[i].pos[4]-vy0;
        vz=(*p)[i].pos[5]-vz0;

        float r;
        if(alpha<0) r = sqrt((x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0)+0.05*0.05*1e-6); // Mpc
        else r = sqrt((x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0));
        //float r = sqrt((*p)[i].r2+1e-10);
        //if((r*1000.-rin) > 0 && (r*1000.-(rin+rout))<0 && ((*p)[i].ke-(*p)[i].pe) < 0){
        if((r*1000.-rin) > 0 && (r*1000.-(rin+rout))<0){
            th=acos((z-cz0)/r);
            phi=atan2((y-cy0),(x-cx0));
            vr = vx*cos(phi)*sin(th)+vy*sin(phi)*sin(th)+vz*cos(th);
            vt = vx*cos(phi)*cos(th)+vy*sin(phi)*cos(th)-vz*sin(th);
            vp = -vx*sin(phi)+vy*cos(phi);
	    sigma0+=(vr-vrbar)*(vr-vrbar)*pow(r*1000.,alpha);
	    //sigma0+=(vt-vtbar)*(vt-vtbar)*pow(r*1000.,alpha);
	    //sigma0+=(vp-vpbar)*(vp-vpbar)*pow(r*1000.,alpha);
	    //sigma0+=(vx-vxbar)*(vx-vxbar)*pow(r*1000.,alpha);
	    //sigma0+=(vy-vybar)*(vy-vybar)*pow(r*1000.,alpha);
	    //sigma0+=(vz-vzbar)*(vz-vzbar)*pow(r*1000.,alpha);
	}
    }
    //sigma0 = sqrt(sigma0/(3.*npart)); 
    sigma0 = sqrt(sigma0/(npart)); 
    //printf("--------------------------------- npart = %f, sigma0 = %f \n",npart,sigma0);

    return sigma0;
}
