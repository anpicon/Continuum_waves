#!/bin/env python3

import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import subprocess
import math

path = os.getcwd()
#----------------------------------------------------------------------------------------#
def read_cubef(fich):
 with open(fich,"r") as f:
  i=int(0)
  Text=[]
  wf=[]
  natoms=int(3)
  for line in f:
   if(i<(natoms+6)):
    Text.append(line)
   s=line.split()
   if(i==2):
    natoms=int(s[0])
   if(i>=(natoms+6)):
    wf.append(s)
   i=i+1
 return Text,wf
 
#----------------------------------------------------------------------------------------#
def write_cubef(fich,Text,wf):
 with open(fich, 'w') as f:
  for item in Text:
   f.write("%s" % item)
  for item in range(0,len(wf)):
   for item2 in wf[item]:
    if (abs(item2) > 1.e-22):
     s  = "%."+str(5)+"E"
     ss = s % (item2)
     f.write("%s " % ss)
   f.write("\n")
 f.close()
 #print(Text)

#----------------------------------------------------------------------------------------#
def read_cubes(folder,list):
 wfs=[]
 for item in list:
  text,cube=read_cubef(folder+item+".cube")
  wfs.append(cube)
 return text,wfs
 
#----------------------------------------------------------------------------------------#
def list2np(wfs):
 n1=len(wfs)
 n2=len(wfs[0])
 n3=len(wfs[0][0])
 phi=np.zeros(n1*n2*n3,dtype=np.float32).reshape(n1,n2,n3)
 for i in range(0,n1):
  for j in range(0,n2):
   for k in range(0,len(wfs[i][j])):
    phi[i][j][k]=float(wfs[i][j][k])
 
 return phi
 
#----------------------------------------------------------------------------------------#
def NORB(diag,root,wincube):#,offd):
 nstate=len(diag)
 nacorb=len(diag[0][:][0])
 
 ndia=np.zeros(nstate*nacorb,dtype=np.float32).reshape(nstate,nacorb)
 
 nactel=[]
 for i in range(0,nstate):
  sum=0.
  for j in range(0,nacorb):
   ndia[i][j]=diag[i,j,j]
   sum=sum+ndia[i][j]
  nactel.append(sum)
 print("# active elec: (per state)")
 print(nactel)
 
 #diagonalizing the density matrix
 occ,u = np.linalg.eig(diag[root-1,:,:])
 print("Root: " + str(root))
 print("Occ:  " + str(occ))
 
 print("Reading wavefunctions...")
 list = np.linspace(wincube[0], wincube[1],wincube[1]-wincube[0]+1).astype(int)
 print("# cube files: ")
 print(list)
 text,wfs=read_cubes(molcasfolder,list.astype(str))

 #constructing and writing the natural orbitals for a specific root
 phi=list2np(wfs)
 n1=len(wfs[0])
 n2=len(wfs[0][0])
 start = int(wincube[1]-nacorb) #we look for the 1st active orbital
 
 occ2 = np.zeros(start,dtype=np.float32)
 for i in range(0,start):
  occ2[i] = 2.
 occ2 = np.append(occ2,occ)
 np.savetxt("occ.txt", occ2, delimiter=',')
 
 for i in range(0,nacorb):
  norb=np.zeros(n1*n2,dtype=np.float32).reshape(n1,n2)
  for j in range(0,len(u[:,i])):
   norb+=u[j,i]*phi[j+start]
  s= "n_orbital_"+str(i+start+1)+".1.cube"
  print("Creating cube file: " + s)
  write_cubef(s,text,norb)
 
 for i in range(0,start):
  s= "n_orbital_"+str(i+1)+".1.cube"
  print("Copying cube file:  " + s)
  write_cubef(s,text,phi[i])
 
 #noff=np.zeros(nstate*nstate*nacorb*nacorb,dtype=np.float32).reshape(nstate*nstate,nacorb,nacorb)
 #count=int(0)
 #for i in range(0,nstate):
 # for j in range(i+1,nstate):
 #  noff[nstate*i+j]=offd[count,:,:]
 #  count=count+1
 
 return occ #ndia,noff
#----------------------------------------------------------------------------------------#

if (len(sys.argv)!=5):
 print("1st arg: molcas folder (h5 file)")
 print("2nd arg: select root")
 print("3nd and 4th arg: range cube files")
 quit()

molcasfolder = sys.argv[1]
iroot        = int(sys.argv[2])
wincube      = np.zeros(2,dtype=np.int)
wincube[0]   = int(sys.argv[3])
wincube[1]   = int(sys.argv[4])

print("Reading occupation numbers...")
with h5py.File(molcasfolder+"input.rasscf.h5","r") as f:
 diag=f['DENSITY_MATRIX'][:]
 #offd=f['TRANSITION_DENSITY_MATRIX'][:]

nstate=len(diag)
nacorb=len(diag[0][:][0])

print("# states     : " + str(nstate))
print("# active orb : " + str(nacorb))

#creating natural orbitals for electron density
occ = NORB(diag,iroot,wincube)


