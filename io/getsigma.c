#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../universal_constants.h"
#include "../config_vars.h"
#include "getsigma.h"

float getsigma(struct potential **p,int64_t num_p,float *cen,float *sigma,float rin,float eps){

    float x,y,z,er,et,ep;
    float vx,vy,vz,vr,vt,vp;
    float th=0,phi=0;
    float vxbar=0,vybar=0,vzbar=0,vrbar=0,vtbar=0,vpbar=0;
    float sigma0=0;

    float npart=0;
    float r2=0;
    float cx0,cy0,cz0;
    float vx0,vy0,vz0;
    float sr0,st0,sp0;
    int64_t i;

    cx0=cen[0];   cy0=cen[1];   cz0=cen[2]; 
    vx0=cen[3];   vy0=cen[4];   vz0=cen[5]; 
    sr0=0; st0=0; sp0=0;

    // particles sort by r2... 
    for(i=0; i<num_p ;i++){
        x=(*p)[i].pos[0];
        y=(*p)[i].pos[1];
        z=(*p)[i].pos[2];
	vx=(*p)[i].pos[3]-vx0;
	vy=(*p)[i].pos[4]-vy0;
	vz=(*p)[i].pos[5]-vz0;

        float r = sqrt((x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0));
        th=acos((z-cz0)/r);
        phi=atan2((y-cy0),(x-cx0));
        vr = vx*cos(phi)*sin(th)+vy*sin(phi)*sin(th)+vz*cos(th);
        vt = vx*cos(phi)*cos(th)+vy*sin(phi)*cos(th)-vz*sin(th);
        vp = -vx*sin(phi)+vy*cos(phi);

        //float r = sqrt((*p)[i].r2+1e-10);
        //if((r*1000.-rin) > 0 && (r*1000.-(rin+eps))<0 && ((*p)[i].ke-(*p)[i].pe) < 0){
        if((r*1000.-rin) > 0 && (r*1000.-(eps))<0){
           // error: need to calculate vx,vy,vz relative to the bulk motion... ... 
           // --------------------------------------------------------------------
           //vr=((x-cx0)*vx+(y-cy0)*vy+(z-cz0)*vz)/r;
           //vrbar+=vr;
           vxbar+=vx; 
           vybar+=vy; 
           vzbar+=vz; 
           vrbar+=vr;
           vtbar+=vt;
           vpbar+=vp;
           npart++;
        }
    }
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

        float r = sqrt((x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0));
        th=acos((z-cz0)/r);
        phi=atan2((y-cy0),(x-cx0));
        vr = vx*cos(phi)*sin(th)+vy*sin(phi)*sin(th)+vz*cos(th);
        vt = vx*cos(phi)*cos(th)+vy*sin(phi)*cos(th)-vz*sin(th);
        vp = -vx*sin(phi)+vy*cos(phi);

        //float r = sqrt((*p)[i].r2+1e-10);
        //if((r*1000.-rin) > 0 && (r*1000.-(rin+eps))<0 && ((*p)[i].ke-(*p)[i].pe) < 0){
        if((r*1000.-rin) > 0 && (r*1000.-(eps))<0){
            //vr=((x-cx0)*vx+(y-cy0)*vy+(z-cz0)*vz)/r;
	    //sigma0+=(vr-vrbar)*(vr-vrbar);
	    // sigma0+=(vx-vxbar)*(vx-vxbar) + (vy-vybar)*(vy-vybar) + (vz-vzbar)*(vz-vzbar);

        /*   ZXY 2022.07.08: sigma0 is actually the dispersion of velocity in central frame of NGC1052 !
                             Because: vx = v[i]_x - vc_x,  so, vxbar = v[i]_x_(ave) - vc_x. 
                             So, vx - vxbar = v[i]_x - v[i]_x_(ave).  This is obviously calculating the dispersion of absolute
                             velocity. But actually, it proves that they differ little ... */
            sigma0 += vx*vx + vy*vy + vz*vz;
            sr0+=(vr-vrbar)*(vr-vrbar);
            st0+=(vt-vtbar)*(vt-vtbar);
            sp0+=(vp-vpbar)*(vp-vpbar);
	}
    }
    sigma0 = sqrt(sigma0/(3.*npart));

    // ZXY 2022.04.20: sigma0 is actually sigma_1D, print it !
    if( (*p)[1].type == 4 && (*p)[2].type == 4 && (*p)[3].type == 4 )
      printf("\n >>>>>>>>> ZXY: sigma_1D_star = %f !!!!!!!!!!!!!!!!!!!!!\n", sigma0);
    if( (*p)[1].type == 1 && (*p)[2].type == 1 && (*p)[3].type == 1 )
      printf("\n\n >>>>>>>>> ZXY: sigma_1D_halo = %f !!!!!!!!!!!!!!!!!!!!!\n\n", sigma0);

    sr0 = sqrt(sr0/npart);
    st0 = sqrt(st0/npart);
    sp0 = sqrt(sp0/npart);

    sigma[0]=sr0;
    sigma[1]=st0;
    sigma[2]=sp0;
    //printf("-------- npart = %f, sigma0 = %f, sr0 = %f, st0 = %f, sp0 = %f, validation %f==%f\n",npart,sigma0,sr0,st0,sp0,sigma0,sqrt((sr0*sr0+st0*st0+sp0*sp0)/3));

    return sigma0;
}
