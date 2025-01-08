import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.interpolate
from scipy.interpolate import interp1d
import matplotlib.patches as mpatches
import math
from scipy.misc import derivative
from math import log
from math import sqrt
import matplotlib.cm as cm
import pickle
import pprint
import os
import sys
import time
import glob
import re
import shutil
import subprocess


folder_name = os.listdir('../../') 
folder_name = [folder for folder in folder_name if 'DM' in folder]
folder_name.sort()
print(folder_name)

j=0

with open("sh-all.sh", 'w') as runall:
  for folder in folder_name:
      runall.write('sh ' + folder + '.sh\n')
      j+=1
      if j%2==0:
        runall.write('sleep 135\n')
      with open(folder + '.sh','w') as sh:
         sh.write('rm ../../'+folder+'/output/boundmass.txt\n')
         sh.write('rm ../../'+folder+'/output/dataProfile/*\n')
         for i in range(60):
            sh.write('nohup srun --mem=10gb --pty ./modRS ../../' +folder+ '/output snapshot {:03} -l &\n'.format(i))
      sh.close()
runall.close




