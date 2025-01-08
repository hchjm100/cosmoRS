#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>
#include <string.h>

void getDMpos_SIDM(){

  char fpath[500];
  char tardir[500];
  char cmd[500];

  char type[5][10];
  int snap[10][10];
  int Orb[10];
  char trash[500];

  FILE *F;

  sprintf(type[0], "SIDM10");
  sprintf(type[1], "SIDM30");
  sprintf(type[2], "SIDM60");

  Orb[0] = 5;  
  Orb[1] = 7;  
  Orb[2] = 3;

  snap[0][0] = 0;   snap[0][1] = 5;    snap[0][2] = 16;   snap[0][3] = 28;  snap[0][4] = 35;  snap[0][5] = 37;  snap[0][6] = 39;
  snap[1][0] = 0;   snap[1][1] = 14;   snap[1][2] = 25;   snap[1][3] = 33;  snap[1][4] = 35;  snap[1][5] = 38;  snap[1][6] = 99;
  snap[2][0] = 0;   snap[2][1] = 10;   snap[2][2] = 22;   snap[2][3] = 35;  snap[2][4] = 42;  snap[2][5] = 44;  snap[2][6] = 46;

  int i, j, k, O, n, ct;


  sprintf(tardir, "/home/storage0/users/zhangxy/CraterII/DMpos-SIDM");
  sprintf(cmd, "rm %s/*", tardir);
  system(cmd);

  for(k=0; k<=2; k++){
    for(n=0; n<=6; n++){

      sprintf(fpath, "/home/storage0/users/zhangxy/CraterII/%s-O%d/output/dataProfile", type[k], Orb[k]);
      sprintf(cmd, "cp %s/DM_positions_%03d.txt  %s/DMpos-%s-O%d-%03d.txt", fpath, snap[k][n], tardir, type[k], Orb[k], snap[k][n]);
      system(cmd);

    }
  }

}



