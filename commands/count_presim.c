#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>
#include <string.h>

void count_presim(){

  char fpath[500];
  char runname[500];
  char cmd[500];
  char type[2][10];
  char trash[500];
  FILE *F;

  sprintf(type[0], "CDM");
  sprintf(type[1], "SIDM");

  int i, j, k, n, ct;
  int ps[2] = {500, 3000};

  
  for(i=1; i<=3; i++)
    {
    for(j=0; j<=1; j++)
      {
      for(k=0; k<=1; k++)
        {
          sprintf(fpath, "/home/storage0/users/zhangxy/CraterII/presim-%s-ps%d-E%d/output/dataProfile/", type[k], ps[j], i);
          sprintf(runname, "presim-%s-ps%d-E%d.ct", type[k], ps[j], i);
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

}



