#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>
#include <string.h>

void getDMprof_CDM(){

  char fpath[500];
  char tardir[500];
  char runname[500];
  char cmd[500];

  char type[2][10];
  int snap[10];
  char trash[500];

  FILE *F;

  sprintf(type[0], "CDM");
  sprintf(type[1], "SIDM");

  snap[0] = 0;
  snap[1] = 7;
  snap[2] = 17;
  snap[3] = 26;
  snap[4] = 35;
  snap[5] = 45;

  int i, j, k, O, n, ct;


  sprintf(tardir, "/home/storage0/users/zhangxy/CraterII/DMprof-CDM");
  sprintf(cmd, "rm %s/*", tardir);
  system(cmd);

  for(k=0; k<1; k++){
   for(O=2; O<=2; O++){
    for(n=0; n<=5; n++){

      sprintf(fpath, "/home/storage0/users/zhangxy/CraterII/%s-ps760-E2-O%d/output/dataProfile", type[k], O);
      sprintf(cmd, "cp %s/DMprof_snap%03d.txt  %s/DMprof-%s-O%d-%03d.txt", fpath, snap[n], tardir, type[k], O, snap[n]);
      system(cmd);

    }
   }
  }

}



