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
    float vxbar[nticks-1];
    float vybar[nticks-1];
    float vzbar[nticks-1];
    float vrbar[nticks-1];
    float vtbar[nticks-1];
    float vpbar[nticks-1];
    int thisnt = nticks;
    float sigx[thisnt-1];
    float sigy[thisnt-1];
    float sigz[thisnt-1];
    float counts[nticks-1];
/* 
    memset(vpfr,  0.0, sizeof(float)*(nticks-1));
    memset(vpft,  0.0, sizeof(float)*(nticks-1));
    memset(vpfp,  0.0, sizeof(float)*(nticks-1));
*/
    memset(vxbar, 0.0, sizeof(float)*(nticks-1));
    memset(vybar, 0.0, sizeof(float)*(nticks-1));
    memset(vzbar, 0.0, sizeof(float)*(nticks-1));
    memset(vrbar, 0.0, sizeof(float)*(nticks-1));
    memset(vtbar, 0.0, sizeof(float)*(nticks-1));
    memset(vpbar, 0.0, sizeof(float)*(nticks-1));
    memset(sigx, 0.0, sizeof(float)*(nticks-1));
    memset(sigy, 0.0, sizeof(float)*(nticks-1));
    memset(sigz, 0.0, sizeof(float)*(nticks-1));
    memset(counts, 0.0, sizeof(float)*(nticks-1));


    int i = 0, ni = 0, ibin =0;
    float r = 0.0;

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
        r = sqrt((x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0));
        th=acos((z-cz0)/r);
        phi=atan2((y-cy0),(x-cx0));
        vr = vx*cos(phi)*sin(th)+vy*sin(phi)*sin(th)+vz*cos(th);
        vt = vx*cos(phi)*cos(th)+vy*sin(phi)*cos(th)-vz*sin(th);
        vp = -vx*sin(phi)+vy*cos(phi);

        r = r*1000; // kpc
        ibin=0;
        while(r>rpf[ibin+1] && ibin<(nticks-1)) ibin++;
        vxbar[ibin]+=vx;
        vybar[ibin]+=vy;
        vzbar[ibin]+=vz;
        vrbar[ibin]+=vr;
        vtbar[ibin]+=vt;
        vpbar[ibin]+=vp;
        counts[ibin]++;
    }
    for(ni=0;ni<nticks-1;ni++){
        if(counts[ni]!=0){
          vxbar[ni]/=counts[ni]; vybar[ni]/=counts[ni]; vzbar[ni]/=counts[ni]; 
          vrbar[ni]/=counts[ni]; vtbar[ni]/=counts[ni]; vpbar[ni]/=counts[ni]; 
        }

    // 2024.06.27 ZXY debug
        if(ni==0)
	  printf(" !~~~~~~~~~~~~~~~~~ZXY debug: count0 = %f !!!!!!!!!!!! \n", counts[ni]);
    }

    // #####################################################
    // Compute velocity dispersions
    for(ni=0; ni<nticks-1; ni++){
	sigz[ni] = 0.0;
	vpfr[ni] = 0.0;
	vpft[ni] = 0.0;
	vpfp[ni] = 0.0;
    }
    for(i=0; i<num_p; i++){
	x=(*p)[i].pos[0];
	y=(*p)[i].pos[1];
	z=(*p)[i].pos[2];
        vx=(*p)[i].pos[3]-vx0;
        vy=(*p)[i].pos[4]-vy0;
        vz=(*p)[i].pos[5]-vz0;

        r = sqrt((x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0));
        th=acos((z-cz0)/r);
        phi=atan2((y-cy0),(x-cx0));
        vr = vx*cos(phi)*sin(th)+vy*sin(phi)*sin(th)+vz*cos(th);
        vt = vx*cos(phi)*cos(th)+vy*sin(phi)*cos(th)-vz*sin(th);
        vp = -vx*sin(phi)+vy*cos(phi);

        r = r*1000; // kpc
        ibin=0;
        while(r>rpf[ibin+1]&&ibin<(nticks-1)) ibin++;
        //vpf[ibin]+=pow(vx-vxbar[ibin],2)+pow(vy-vybar[ibin],2)+pow(vz-vzbar[ibin],2);
/*
	vpfr[ibin]+=pow(vr-vrbar[ibin],2);
        vpft[ibin]+=pow(vt-vtbar[ibin],2);
        vpfp[ibin]+=pow(vp-vpbar[ibin],2);
*/
	sigx[ibin]+=vx*vx;
	sigy[ibin]+=vy*vy;
	sigz[ibin]+=pow(vz,2);

	if(ibin==0){
	  printf("~~~~~~~~~~ZXY debug: ibin vpfr = %f, vpft = %f, vpfp = %f ~~~~~~~~~\n", sigx[ibin], sigy[ibin],sigz[ibin]);
	  sigx[ibin] = 0; sigy[ibin] = 0;  sigz[ibin] = 0;
	}
    }

    printf("~~~~~~~~~~ZXY debug: vpfr = %f, vpft = %f, vpfp = %f ~~~~~~~~~\n", sigx[0], sigy[0], sigz[0]);

    ibin = 0;
    int nin = 0;
    for(nin=0; nin<nticks-1; nin++){ 
        if(counts[nin] > 20){ 
            //vpf[ni] = sqrt(vpf[ni]/(3.*counts[ni]));
            vpfr[nin] = sqrt(sigx[nin]/counts[nin]); 
            vpft[nin] = sqrt(sigy[nin]/counts[nin]); 
            vpfp[nin] = sqrt(sigz[nin]/counts[nin]); 
        }
        else{
            vpfr[nin] = 0;
            vpft[nin] = 0;
            vpfp[nin] = 0;
        }
    }

    return 0;
}
