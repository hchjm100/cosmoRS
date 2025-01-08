#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>


static double All_G = 4.30091727e-6;


double phi_disk_ZXY(double M_d, double a_d, double b_d, double R, double z){

  return -All_G * M_d / sqrt( R*R + pow( a_d + sqrt(b_d*b_d + z*z), 2) );

}

double phi_Hern_ZXY(double M_h, double r_h, double r){

  return -All_G * M_h / (r_h + r);

}


// 2023.04.20 ZXY: Actually all the forces do not diverge as r -> 0, therefore no need of softening ......
//                 I just take the limit of r -> 0 ......
double phi_NFW_ZXY(double rhos, double rs, double r){

    return -4*M_PI * All_G * pow(rs,3) * rhos * log(1+r/rs) / r;

}


  // DY: external potential, following Volker's suggestion 
  // Pos: comoving coordinates in units h^-1 kpc
  // physical pos: Pos*a/h = Pos*All.Time/All.HubbleParam
  // printf("\n All.Time %g h: %g \n",All.Time,All.HubbleParam);


  // 2023.03.31 ZXY: Begin MW potential ......
  // 2024.03.12 ZXY: Add a switch of external potential ......
  // 2024.04.17 ZXY: Revise MW potential according to McMillan (2017) ...
  // 2024.04.18 ZXY: According to McMillan (2017), there are other two gas discs, named d3 and d4 here ...


static  double M_d1 = 3.518583772020568e10,   a_d1 = 2.50,   b_d1 = 0.3;
static  double M_d2 = 1.048684487943493e10,   a_d2 = 3.02,   b_d2 = 0.9;
static  double M_d3 = 1.1e10,      a_d3 = 7,      b_d3 = 0.085;
static  double M_d4 = 0.12e10,     a_d4 = 1.5,    b_d4 = 0.045;

static  double M_Hern = 0.923e10,   r_Hern = 1.3;
static  double rhos_NFW = 8.54e6,  rs_NFW = 19.6;

double phi_MW(double x, double y, double z){

  double rad_ZXY = sqrt( x*x + y*y + z*z );
  double Rad_ZXY = sqrt( x*x + y*y );

  double phi_d1 = phi_disk_ZXY(M_d1, a_d1, b_d1, Rad_ZXY, z);
  double phi_d2 = phi_disk_ZXY(M_d2, a_d2, b_d2, Rad_ZXY, z);
  double phi_d3 = phi_disk_ZXY(M_d3, a_d3, b_d3, Rad_ZXY, z);
  double phi_d4 = phi_disk_ZXY(M_d4, a_d4, b_d4, Rad_ZXY, z);

  double phi_H  = phi_Hern_ZXY(M_Hern, r_Hern, rad_ZXY);
  double phi_N  = phi_NFW_ZXY(rhos_NFW, rs_NFW, rad_ZXY);

  return (phi_d1 + phi_d2 + phi_d3 + phi_d4 + phi_H + phi_N);

}

