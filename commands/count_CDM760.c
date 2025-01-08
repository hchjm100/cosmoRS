#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>
#include <string.h>

void count_CDM760(){

  char fpath[500];
  char runname[500];
  char cmd[500];
  char type[2][10];
  char trash[500];
  FILE *F;

  sprintf(type[0], "CDM");
//  sprintf(type[1], "SIDM");

  int i, j, k, l, n, ct;
  int ps[1] = {760};

  
  for(i=2; i<=2; i++)
    {
    for(j=0; j<1; j++)
      {
      for(k=0; k<1; k++)
        {
        for(l=2; l<=2; l++)
          {
          sprintf(fpath, "/home/storage0/users/zhangxy/CraterII/%s-ps%d-E%d-O%d/output/dataProfile/", type[k], ps[j], i, l);
          sprintf(runname, "%s-ps%d-E%d-O%d.ct", type[k], ps[j], i, l);
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

}



