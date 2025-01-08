#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../universal_constants.h"
#include "../config_vars.h"
#include "getMprof.h"

float getMprof(struct potential **p,int64_t num_p,float *cen,float *mr,float *rpf,const int nticks){

    double  tor = (double)SCALE_NOW/h0;
    double  tom = 1.0/h0;
    float x,y,z;
    float vx,vy,vz;
    float vx1,vy1,vz1;

    float npart=0;
    float cx0,cy0,cz0;
    float cx1,cy1,cz1;
    float rden,mr0;
    float sw[nticks-1]; 
    //memset(sw,  0, sizeof(float)*(nticks-1)); 
    memset(mr, 0, sizeof(float)*(nticks-1)); 
    int64_t i,j;

    cx0=cen[0]; cy0=cen[1]; cz0=cen[2]; 
    // Mass unit: solar mass!
    //for(int ni=1;ni<nticks;ni++){
    //}
    float r=0;
    for(i=0; i<num_p ;i++){
        x=(*p)[i].pos[0];
        y=(*p)[i].pos[1];
        z=(*p)[i].pos[2];
        r = sqrt((x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0));
        r = r*1000; // kpc
        // float rlow=rpf[0];
        // float rhigh=rpf[1];
        int ibin=0;
        while(r>rpf[ibin+1] && ibin<(nticks-1)) ibin++;
        //mr[ibin]+=sw[ibin];
	// 2024.04.18 ZXY: Support different particle mass ......
	if(ibin<nticks-1){
          for(j=ibin;j<nticks;j++) mr[j]+= PART_MASS_ZXY[(*p)[i].type];
	}
    }

    return 0;
}
