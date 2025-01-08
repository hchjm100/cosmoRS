#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

extern double All_G;

extern double M_d1, a_d1, b_d1;
extern double M_d2, a_d2, b_d2;
extern double M_d3, a_d3, b_d3;
extern double M_d4, a_d4, b_d4;
extern double M_Hern, r_Hern;
extern double rhos_NFW, rs_NFW;

double phi_MW(double x, double y, double z);
double phi_disk_ZXY(double M, double a, double b, double Rad_ZXY, double z);
double phi_Hern_ZXY(double M_Hern, double r_Hern, double rad_ZXY);
double phi_NFW_ZXY(double rhos_NFW, double rs_NFW, double rad_ZXY);
