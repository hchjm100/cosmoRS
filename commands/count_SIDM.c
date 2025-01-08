#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>
#include <string.h>

void count_SIDM(){

  char fpath[500];
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

  
      for(k=1; k<=4; k++)
        {
        for(l=2; l<=11; l++)
          {
          sprintf(fpath, "/home/storage0/users/zhangxy/CraterII/%s-O%d/output/dataProfile/", type[k], l);
          sprintf(runname, "%s-O%d.ct", type[k], l);
          sprintf(cmd, "ls %s |grep txt > %s", fpath, runname);
          system(cmd);

          if(!(F = fopen(runname,"r")))
            printf("Opening %s failed !!!", runname);


          for(ct=0; !feof(F) && fgets(trash,499,F); ct++);


          fclose(F);
          printf("%s: count = %d ... \n", runname, ct);


          sprintf(cmd, "rm %s", runname);
          system(cmd);
          }
        }

}



