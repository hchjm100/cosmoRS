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


folder_name = glob.glob('*DM*')
folder_name = sorted(folder_name)
print(folder_name)

subprocess.run("rm -r dataExtract/*", shell=True)

for folder in folder_name:
    subprocess.run(['mkdir', './dataExtract/' +folder])
    subprocess.run('cp ./' +folder + '/output/dataProfile/DMprof* ./dataExtract/' +folder+ '/', shell=True)
    subprocess.run('cp ./' +folder + '/output/dataProfile/prof_star_snapshot* ./dataExtract/' +folder+ '/', shell=True)
    subprocess.run('cp ./' +folder + '/output/boundmass.txt ./dataExtract/' +folder+ '/', shell=True)


