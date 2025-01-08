#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/stat.h>
#include <sys/types.h>

void mvfile(){

  char fpath[500];
  char runname[500];
  char type[10][10];

  sprintf(type[0], "CDM");
  sprintf(type[1], "SIDM10");
  sprintf(type[2], "SIDM30");
  sprintf(type[3], "SIDM100");

  ofstream in;

  int i, j, k, l, n;

  int ps[2] = {760,3000};

            in.open("mvsnap.sh");

            for(n=1; n<=55; n++)
              {
                in << Form("mv snapshot_%03d snapshot_%03d", n+22, n) << endl;
              }

            in.close();

}


/*
  for(i=1; i<=3; i++)
    {
      sprintf(fpath, "/home/storage0/users/zhangxy/CraterII/presim-CDM-ps%d-E%d/output", 500, i);
      sprintf(runname, "presim-SIDM-ps%d-E%d.sh", ps1, i);
      in.open(runname);
      in << Form("rm ") << fpath << Form("/boundmass.txt") << endl;
      for(n=0; n<=25; n++)
      {
        in << Form("nohup ./modRS ") << fpath << Form(" snapshot %03d &", n) << endl;
      }
      in.close();
    }

  for(i=1; i<=3; i++)
    {
      sprintf(fpath, "/home/storage0/users/zhangxy/CraterII/presim-CDM-ps%d-E%d/output", 500, i);
      sprintf(runname, "presim-CDM-ps%d-E%d.sh", ps2, i);
      in.open(runname);
      in << Form("rm ") << fpath << Form("/boundmass.txt") << endl;
      for(n=0; n<=25; n++)
      {
        in << Form("nohup ./modRS ") << fpath << Form(" snapshot %03d &", n) << endl;
      }
      in.close();
    }

  for(i=1; i<=3; i++)
    {
      sprintf(fpath, "/home/storage0/users/zhangxy/CraterII/presim-CDM-ps%d-E%d/output", 500, i);
      sprintf(runname, "presim-SIDM-ps%d-E%d.sh", ps2, i);
      in.open(runname);
      in << Form("rm ") << fpath << Form("/boundmass.txt") << endl;
      for(n=0; n<=25; n++)
      {
        in << Form("nohup ./modRS ") << fpath << Form(" snapshot %03d &", n) << endl;
      }
      in.close();
    }
*/

