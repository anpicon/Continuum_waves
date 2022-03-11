#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import subprocess
from subprocess import Popen, PIPE
import shutil

path = os.getcwd()
program = "/home/jgv/bin/grid2cube.exe"

#----------------------------------------------------------------------------------------#
 
if (len(sys.argv)!=4):
 print("1st arg: lus file")
 print("2st arg: label output file")
 print("3st arg: # output files")
 quit()

infile = sys.argv[1]
fout   = sys.argv[2]
nmo    = int(sys.argv[3])

for i in range(0,nmo):
 #outputfile = fout + "_" + str(i+1) + ".cube"
 outputfile = str(i+1) + ".cube"
 proc = subprocess.Popen([program,infile,outputfile],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
 proc.communicate(input=str(i+1).encode())
 proc.wait()

