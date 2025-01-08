#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../universal_constants.h"
#include "../config_vars.h"
#include "getprof.h"

float getprof(struct potential **p,int64_t num_p,float *cen,float *rho,float *rpf,const int nticks){

    double  tor = (double)SCALE_NOW/h0;

    float x,y,z;
    float vx,vy,vz;
    float vx1,vy1,vz1;

    float npart=0;
    float cx0,cy0,cz0;
    float cx1,cy1,cz1;
    float rden,rho0;
    float sw[nticks-1]; 
    memset(sw,  0, sizeof(float)*(nticks-1)); 
    memset(rho, 0, sizeof(float)*(nticks-1)); 
    int64_t i, j, k, l, m, n;

    // ZXY 2022.04.18: Add a process of calculating 3D and 2D projected star density ......
/*
    if (num_p < 10000)      // For stars only ......
    {  

      double size = 20;
      int nmesh = 100;
      double d;
      d = 2 * size / ((double)(nmesh));

      printf("\n >>>>>>> ZXY: Constructing 3D star density of (%f kpc)^3, dx = %f kpc\n", size, d);

      double rho2xy[nmesh+1][nmesh+1], rho2xz[nmesh+1][nmesh+1], rho2yz[nmesh+1][nmesh+1];
      double rho3[nmesh+1][nmesh+1][nmesh+1];
      double axis[nmesh+1];

      for( i=0; i<=nmesh; i++ )
        axis[i] = -size + d/2 + ((double)(i))*d;

      printf("\n Test mesh = %f, %f, %f, %f, %f .....\n", axis[0], axis[10], axis[30], axis[50], axis[100]);

      for( i=0; i<=nmesh; i++ ){
        for( j=0; j<=nmesh; j++ ){
            axis2[i][j] = 0;
          for( k=0; k<=nmesh; k++ ){
            axis3[i][j][k] = 0;
          }
        }
      }


    }
*/

    cx0=cen[0]; cy0=cen[1]; cz0=cen[2]; 
    // Mass unit: solar mass!
    for(int ni=1;ni<nticks;ni++){
       // 2024.04.18 ZXY: Support different particle mass ......
       sw[ni-1]=PART_MASS_ZXY[(*p)[0].type]/(4./3.*M_PI*(pow(rpf[ni],3)-pow(rpf[ni-1],3)));
    }

    float r=0;
    for(i=0; i<num_p ;i++){

       // ZXY 2022.02.12: Exclude unbounded particles ?
       if( (*p)[i].pe < (*p)[i].ke) continue;

        x=(*p)[i].pos[0];
        y=(*p)[i].pos[1];
        z=(*p)[i].pos[2];
        r = 1000*sqrt((x-cx0)*(x-cx0)+(y-cy0)*(y-cy0)+(z-cz0)*(z-cz0));
       // float rlow=rpf[0];
       // float rhigh=rpf[1];
        int ibin=0;
        while(r>rpf[ibin+1]&&ibin<(nticks-1)) ibin++;
	if(ibin<nticks-1)
          rho[ibin]+=sw[ibin]; 
    }

    return 0;
}
