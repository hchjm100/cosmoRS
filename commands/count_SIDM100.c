#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>
#include <string.h>

void count_SIDM100(){

  char fpath[500];
  char runname[500];
  char cmd[500];
  char type[10][10];
  char trash[500];
  FILE *F;

  char SIDM100O[10][10];

  sprintf(SIDM100O[0], "7.1");
  sprintf(SIDM100O[1], "7.2");
  sprintf(SIDM100O[2], "7.3");
  sprintf(SIDM100O[3], "8");


  int i, j, k, l, n, ct;

  
        for(l=0; l<=3; l++)
          {
          sprintf(fpath, "/home/storage0/users/zhangxy/CraterII/SIDM100-O%s/output/dataProfile/", SIDM100O[l]);
          sprintf(runname, "SIDM100-O%s.ct", SIDM100O[l]);
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



