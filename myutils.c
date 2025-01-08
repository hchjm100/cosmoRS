#include <math.h>
#include "universal_constants.h"
#include "myutils.h"

double arrMOverR[100000];

double getMOR(int64_t i){
   // lookup table: 10 pc resolution
   // 800 kpc=800,00 * 10 pc

   double  r  = ((float)i)/(1e5); // 10 pc to Mpc
   double  rhos=0.000160686e19;// e10/e-9
   double  rs = 0.080; // Mpc

   if(r>0.000001) return (4*3.14159265358979*rhos*pow(rs,3)*(-(r/(r + rs)) - log(rs/(r + rs))))/r;
   else return 0;

};

void setMOR(float *cen){
   int64_t i,j;
   float rh=sqrt(cen[0]*cen[0]+cen[1]*cen[1]+cen[2]*cen[2]);//in Mpc
   j=(int)(rh*1e5);
   //memset(arrMOverR, 0, sizeof(arrMOverR)/sizeof(float));
   //arrMOverR=malloc(100000 * sizeof(float));
   for(i=0;i<(100000);i++){ 
       arrMOverR[i]=getMOR(i)-getMOR(j);
       //if(i%10000==0) printf("-----%d,  %f\n",i,arrMOverR[i]);
   }
};

// Tidal field of the main halo
// Should avoid using this for many particles
double getDUT(int dct,float *pos){

   float  rhos=0.000160686e19;// e10/e-9
   float  rs = 0.080; // Mpc
   double  x,y,z,dutx,duty,dutz;

   x=pos[0]; y=pos[1]; z=pos[2];

   if(fabs(x)>0.000001&&fabs(y)>0.000001&&fabs(z)>0.000001){
      if(dct==0){ dutx= (-4*Gc*3.14159265358979*rhos*pow(rs,3)*x*(2*pow(x,2) + 2*pow(y,2) + 2*pow(z,2) + rs*sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + (pow(rs,2) + pow(x,2) + pow(y,2) + pow(z,2) + 2*rs*sqrt(pow(x,2) + pow(y,2) + pow(z,2)))*log(rs) - (pow(rs,2) + pow(x,2) + pow(y,2) + pow(z,2) + 2*rs*sqrt(pow(x,2) + pow(y,2) + pow(z,2)))*log(rs + sqrt(pow(x,2) + pow(y,2) + pow(z,2)))))/ (pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*pow(rs + sqrt(pow(x,2) + pow(y,2) + pow(z,2)),2)); return dutx; }

      if(dct==1){ duty= (-4*Gc*3.14159265358979*rhos*pow(rs,3)*y*(2*pow(x,2) + 2*pow(y,2) + 2*pow(z,2) + rs*sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + (pow(rs,2) + pow(x,2) + pow(y,2) + pow(z,2) + 2*rs*sqrt(pow(x,2) + pow(y,2) + pow(z,2)))*log(rs) - (pow(rs,2) + pow(x,2) + pow(y,2) + pow(z,2) + 2*rs*sqrt(pow(x,2) + pow(y,2) + pow(z,2)))*log(rs + sqrt(pow(x,2) + pow(y,2) + pow(z,2)))))/(pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*pow(rs + sqrt(pow(x,2) + pow(y,2) + pow(z,2)),2)); return duty; }

      if(dct==2){ dutz= (-4*Gc*3.14159265358979*rhos*pow(rs,3)*z*(2*pow(x,2) + 2*pow(y,2) + 2*pow(z,2) + rs*sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + (pow(rs,2) + pow(x,2) + pow(y,2) + pow(z,2) + 2*rs*sqrt(pow(x,2) + pow(y,2) + pow(z,2)))*log(rs) - (pow(rs,2) + pow(x,2) + pow(y,2) + pow(z,2) + 2*rs*sqrt(pow(x,2) + pow(y,2) + pow(z,2)))*log(rs + sqrt(pow(x,2) + pow(y,2) + pow(z,2)))))/(pow(pow(x,2) + pow(y,2) + pow(z,2),1.5)*pow(rs + sqrt(pow(x,2) + pow(y,2) + pow(z,2)),2)); return dutz; }
      return 0;
   }
   else return 0;
}

