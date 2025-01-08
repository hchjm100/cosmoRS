#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>
#include <string.h>

void getbdmass_CDM_SIDM(){

  char fpath[500];
  char tardir[500];
  char runname[500];
  char cmd[500];
  char type[10][10];
  char trash[500];
  FILE *F;

  sprintf(type[0], "CDM");
  sprintf(type[1], "SIDM10");
  sprintf(type[2], "SIDM30");
  sprintf(type[3], "SIDM50");
  sprintf(type[4], "SIDM100");
  sprintf(type[5], "SIDM70");

  int i, j, k, l, n, ct;


  sprintf(tardir, "/home/storage0/users/zhangxy/CraterII/boundmass-collect");
  sprintf(cmd, "rm %s/*", tardir);
  system(cmd);

      for(k=1; k<=5; k++)
        {
          for(l=2; l<=11; l++)
          {
            sprintf(fpath, "/home/storage0/users/zhangxy/CraterII/%s-O%d/output", type[k], l);
            sprintf(cmd, "cp %s/boundmass.txt  %s/boundmass-%s-O%d.txt", fpath, tardir, type[k], l);
            system(cmd);
          }
        }

}



