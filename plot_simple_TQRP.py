#!/usr/bin/python2.7

'''
minimal plotting tool for pMicro
'''

import sys,math,os,subprocess
#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

varfilename="varMicro.txt"
parfilename="parMicro.txt"
Time=[];health=[];cellvol=[];grate=[];
R=[];R1=[];Q=[];T=[];E=[]
mrp=[];mq=[] ;mt=[] ;me=[] ;
s=[] ;c=[] ;a=[];
tarvol=[]; p_lifespan=[]; 
cell_id=[]; bestguy=[]; ancestors=[];
abund_fittest=[];ndiffgen=[];
av_trm=[]

maxTime=9999999999.
if len(sys.argv)>1:
  if sys.argv[1]=='-maxtime':
    maxTime=float(sys.argv[2])
  
files = [f for f in os.listdir('./')]
title=[f for f in files if 'pMicro0.' in f][0]

with open(varfilename,"r") as fin:
      for line in fin:
        line=line.split()
        
        time = float(line[0])
        if time>maxTime: break
        #if time>500: break
        Time.append( time )
        health.append( float(line[1]) )
        cellvol.append( float(line[2]) )
        R.append(float(line[3]))
        R1.append(float(line[4]))
        Q.append(float(line[5]))
        T.append(float(line[6]))
        mq.append(float(line[7]))
        mt.append(float(line[8]))
        mrp.append(float(line[9]))
        s.append(float(line[10]))
        a.append(float(line[11]))
        grate.append(float(line[12]))
        tarvol.append( float(line[13]) )
        p_lifespan.append( float(line[14]) )
        cell_id.append( int(line[15]) )
        bestguy.append( int(line[16]) )
        ancestors.append( int(line[17]) )
        abund_fittest.append( int(line[18]) )
        ndiffgen.append( int(line[19]) )
        av_trm.append( float(line[20]) )
        
kt=[];at=[];bt=[];ke=[];ae=[];be=[]; 

ktrt0=[];ktrt1=[];ktrt2=[];ktrt3=[];
ktre0=[];ktre1=[];ktre2=[];ktre3=[];
ktrq0=[];ktrq1=[];ktrq2=[];ktrq3=[];
ktrrp0=[];ktrrp1=[];ktrrp2=[];ktrrp3=[];
ktrrr0=[];ktrrr1=[];ktrrr2=[];ktrrr3=[];
geneT=[];geneE=[];geneQ=[];geneRp=[];geneRr=[];
dgeneT=[];dgeneE=[];dgeneQ=[];dgeneRp=[];dgeneRr=[];
Time1=[]
mutinsertion=[];mutdeletion=[];mutdeleterious=[];

with open(parfilename,"r") as fin:
      for line in fin:
        line=line.split()
        
        time = float(line[0])
        if time>maxTime: break
        
        Time1.append(time)
        #if time<200: continue
        #if time>500: break
                
        kt.append( float(line[1]) )
        at.append( float(line[2]) )
        bt.append(float(line[3]))
        
        ktrt0.append(float(line[4]))
        ktrt1.append(float(line[5]))
        
        ktrq0.append(float(line[6]))
        ktrq1.append(float(line[7]))
        
        ktrrr0.append(float(line[8]))
        ktrrr1.append(float(line[9]))
        
        ktrrp0.append(float(line[10]))
        ktrrp1.append(float(line[11]))
                
        geneT.append(int(line[12]))
        geneQ.append(int(line[13]))
        geneRr.append(int(line[14]))
        geneRp.append(int(line[15]))
        
        dgeneT.append(int(line[16]))
        dgeneQ.append(int(line[17]))
        dgeneRr.append(int(line[18]))
        dgeneRp.append(int(line[19]))
        
        mutinsertion.append(float(line[20]))
        mutdeletion.append(float(line[21]))
        mutdeleterious.append(float(line[22]))

######################################
###### ----                ---- ######
###### ---     PLOTTING     --- ######
###### -----              ----- ######
######################################
'''
namefouts="environment.txt"
print "Done reading, now processing..."

if os.path.isfile(namefouts): 
  bla = subprocess.Popen(['tail', '-n', '1', namefouts], stdout=subprocess.PIPE)
  bla = bla.communicate()[0]
  if len(bla)>0:
    initT=float(bla.split(' ')[0])
    indexinitT=Time.index(initT)
    Time_toprint=Time[indexinitT+1:]
    s_toprint=s[indexinitT+1:]
  else:
    Time_toprint=Time
    s_toprint=s
else:
  Time_toprint=Time
  s_toprint=s

tbla=zip(Time_toprint,s_toprint)
with open(namefouts,"a") as fouts:
  for bla in tbla:
    fouts.write(str(bla[0])+' '+str(bla[1])+'\n')
'''

fontsize="x-small"

nsubplots=10
# Two subplots, unpack the axes array immediately
fig, ax = plt.subplots(nsubplots, 1, sharex=True)

subplotnum=0
ax[subplotnum].plot(Time, health, label="health")
ax[subplotnum].plot(Time,R, label="Rp")
ax[subplotnum].plot(Time,R1, label="Rr")
ax[subplotnum].plot(Time,Q, label="Q")
ax[subplotnum].plot(Time,T, label="T")


subplotnum=1
ax[subplotnum].plot(Time, mq, label="mq")
ax[subplotnum].plot(Time, mt, label="mt")
ax[subplotnum].plot(Time, mrp, label="mrp")

subplotnum=2
ax[subplotnum].plot(Time,s,label="S")
ax[subplotnum].plot(Time,a,label="A")

subplotnum=3
#ax[subplotnum].plot(Time1,kt,label="kt")
#ax[subplotnum].plot(Time1,at,label="at")
#ax[subplotnum].plot(Time1,bt,label="bt")
ax[subplotnum].plot(Time1,geneT,label="geneT")
ax[subplotnum].plot(Time1,geneQ,label="geneQ")
ax[subplotnum].plot(Time1,geneRr,label="geneRr")
ax[subplotnum].plot(Time1,geneRp,label="geneRp")
subplotnum=4
ax[subplotnum].plot(Time1,dgeneT,label="dgeneT")
ax[subplotnum].plot(Time1,dgeneQ,label="dgeneQ")
ax[subplotnum].plot(Time1,dgeneRr,label="dgeneRr")
ax[subplotnum].plot(Time1,dgeneRp,label="dgeneRp")
subplotnum=5
ax[subplotnum].plot(Time,grate,label="grate")

subplotnum=6
ax[subplotnum].plot(Time1,ktrt0,label="ktrt0")
ax[subplotnum].plot(Time1,ktrt1,label="ktrt1")
ax[subplotnum].plot(Time1,ktrq0,label="ktrq0")
ax[subplotnum].plot(Time1,ktrq1,label="ktrq1")
subplotnum=7
ax[subplotnum].plot(Time1,ktrrr0,label="ktrrr0")
ax[subplotnum].plot(Time1,ktrrr1,label="ktrrr1")
ax[subplotnum].plot(Time1,ktrrp0,label="ktrrp0")
ax[subplotnum].plot(Time1,ktrrp1,label="ktrrp1")

subplotnum=8
ax[subplotnum].plot(Time1,mutinsertion,label="insert")
ax[subplotnum].plot(Time1,mutdeletion,label="delete")
ax[subplotnum].plot(Time1,mutdeleterious,label="bad")
ax[subplotnum].plot(Time,av_trm,label="av trm")

'''
# Shrink current axis by 20%
box = ax[9].get_position()
ax[9].set_position([box.x0, box.y0, box.width * 0.85, box.height])
# Put a legend to the right of the current axis
ax[9].legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fontsize)
'''
subplotnum=9
s=0.5
ax[subplotnum].plot(Time,cellvol, label="vol")
#ax[subplotnum].plot(Time,tarvol, label="tvol")
ax[subplotnum].plot(Time,p_lifespan, label="life sp.")
#ax[subplotnum].plot(Time,cell_id, label="id bg")
#ax[subplotnum].plot(Time,bestguy, label="pos bg")
#ax[subplotnum].plot(Time,ancestors, label="n anc")
#ax[subplotnum].plot(Time,abund_fittest,label="n fit")
#ax[subplotnum].plot(Time,ndiffgen,label="n diff G")

for n in range(nsubplots):
  # Shrink current axis by 20%
  box = ax[n].get_position()
  ax[n].set_position([box.x0, box.y0, box.width * 0.85, box.height])
  # Put a legend to the right of the current axis
  ax[n].legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fontsize,ncol=2)

plt.suptitle(title,fontsize=12)

plt.savefig("time_plot.png")
plt.savefig("time_plot.svg")
plt.show()  










