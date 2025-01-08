#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../universal_constants.h"
#include "../config_vars.h"
#include "getsigmaprof.h"

float getsigmaprof(struct potential **p,int64_t num_p,float *cen,float *vpfr,float *vpft,float *vpfp,float *rpf,const int nticks){

    float x,y,z,er,et,ep;
    float vx1,vy1,vz1,vr1,vt1,vp1;
    float vx,vy,vz,vr,vt,vp;
    float vx0,vy0,vz0;
    float th=0,phi=0;

    float npart=0;
    float r2=0;
    float cx0,cy0,cz0;
    float cx1,cy1,cz1;
    float rden,v0;
    float sw[nticks-1];
    float vxbar[nticks-1];
    float vybar[nticks-1];
    float vzbar[nticks-1];
    float vrbar[nticks-1];
    float vtbar[nticks-1];
    float vpbar[nticks-1];
    float counts[nticks-1];
    memset(sw,   0, sizeof(float)*(nticks-1));
    memset(vpfr,  0, sizeof(float)*(nticks-1));
    memset(vpft,  0, sizeof(float)*(nticks-1));
    memset(vpfp,  0, sizeof(float)*(nticks-1));
    memset(vxbar, 0, sizeof(float)*(nticks-1));
    memset(vybar, 0, sizeof(float)*(nticks-1));
    memset(vzbar, 0, sizeof(float)*(nticks-1));
    memset(vrbar, 0, sizeof(float)*(nticks-1));
    memset(vtbar, 0, sizeof(float)*(nticks-1));
    memset(vpbar, 0, sizeof(float)*(nticks-1));
    memset(counts,0,sizeof(float)*(nticks-1));
    int64_t i,ni;

    cx0=cen[0];   cy0=cen[1];   cz0=cen[2]; 
    vx0=cen[3];   vy0=cen[4];   vz0=cen[5];

    // #####################################################
    // Compute averaged velocities

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

        r = r*1000; // kpc
        int ibin=0;
        while(r>rpf[ibin+1] && ibin<(nticks-1)) ibin++;
	if(ibin<nticks-1){
          vxbar[ibin]+=vx;
          vybar[ibin]+=vy;
          vzbar[ibin]+=vz;
          vrbar[ibin]+=vr;
          vtbar[ibin]+=vt;
          vpbar[ibin]+=vp;
          counts[ibin]++;
	}
    }
    for(ni=0;ni<nticks-1;ni++){
        if(counts[ni]!=0){
          vxbar[ni]/=counts[ni]; vybar[ni]/=counts[ni]; vzbar[ni]/=counts[ni]; 
          vrbar[ni]/=counts[ni]; vtbar[ni]/=counts[ni]; vpbar[ni]/=counts[ni]; 
        }
    }

    // #####################################################
    // Compute velocity dispersions

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

        r = r*1000; // kpc
        int ibin=0;
        while(r>rpf[ibin+1]&&ibin<(nticks-1)) ibin++;
        //vpf[ibin]+=pow(vx-vxbar[ibin],2)+pow(vy-vybar[ibin],2)+pow(vz-vzbar[ibin],2);
	if(ibin < nticks-1){ 
	  // 2024.06.28 ZXY: Fixed a bug that may cause strange values of dispersion. If no ibin < nticks-1 condition, 
	  // vpfr may leap outside the iteration. That is because the particles outside the last rpf bin will be added 
	  // to vpfr[nticks-1], an element which does not exist, thus cause error on the value of vpfr[0] ......
          vpfr[ibin]+=pow(vr,2);
          vpft[ibin]+=pow(vt,2);
          vpfp[ibin]+=pow(vp,2);
	}
    }
    for(ni=0;ni<nticks-1;ni++){ 
        if(counts[ni] > 0){
            //vpf[ni] = sqrt(vpf[ni]/(3.*counts[ni])); 
            vpfr[ni] = sqrt(vpfr[ni]/counts[ni]); 
            vpft[ni] = sqrt(vpft[ni]/counts[ni]); 
            vpfp[ni] = sqrt(vpfp[ni]/counts[ni]); 
        }
        else{
            vpfr[ni] = 0;
            vpft[ni] = 0;
            vpfp[ni] = 0;
        }
    }

    return 0;
}
