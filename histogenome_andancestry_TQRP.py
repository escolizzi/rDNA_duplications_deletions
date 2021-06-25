#!/usr/bin/python2.7

'''
this script scans the backup files, which contain dumps of the entire population
and prints population stuff, like histgrams of the genomesizes in the population 
and ancestor tracings...
'''

import sys,math,os,string,subprocess
#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

lanc=[]
lancgen=[]
lancT=[];lancQ=[];lancR=[];lancB=[];lancP=[]
lanct=[];lancq=[];lancr=[];lancp=[]

# COLOR MAP #
colorT=(56./255.,146./255.,56./255.)# green (medium blue)
colorQ=(0/255.,109./255.,240/255.)	# blue (lotsa blue)
colorRr=(166/255.,0/255.,0/255.)	# red (no blue)
colorRp=(26/255.,26/255.,30/255.)	# grey - antracitis

color1=(6./255.,6./255.,98/255.)
color2=(128/255.,0/255.,0/255.)
color3=colorT
color4=colorRp

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap



#Transcription - returns only fraction of transcription relative to max multiplied by the fraction of active genes
def Expression_Stuff( lk0,lA,lk1,lhealth, par_max_expression_per_gene, lancAct , lancInact  ):
  fr_=[]
  tr_=[]
  factinact_=[]
  
  req_=[ (x[0]+x[1]*x[2])*x[3] for x in zip(lk0,lA,lk1,lhealth) ]	# requested amount of Rr = (k0+k1*A)*health
  max_=[  par_max_expression_per_gene * (x[0]+x[1]) for x in zip( lancAct,lancInact )  ] # max amount of Rr = par_max_expression_per_gene
  
  for x in zip(req_,max_, lancAct,lancInact):
    
    #if we have zero active genes, then there is a bug somewhere, 
    # because it is impossible that this guy makes in the line of descent
    if x[2]==0.: 
      print "Warning, zero active genes, this is bad!"
    fr_act_inact = x[2]/float(x[2]+x[3])  #fraction of active genes over total act+inact
    factinact_.append( fr_act_inact )
    
    if x[0] <= 0.:
      fr_.append(0.)
      tr_.append(0.)
      
    else:
      ratio = x[0]/x[1] #ratio req_ / max_
      
      
      
      if ratio>1.: 
        fr_.append(1.)
        tr_.append(x[1]* fr_act_inact )
          
      else: 
        fr_.append(ratio)
        tr_.append( x[0]* fr_act_inact )
    
  return fr_,factinact_,tr_


#################################
##### --- PLOTTING DATA --- #####
#################################
def PlotStuff(*arg):
  ltime=arg[0]
  lhist1=arg[1]
  lanc=arg[2]
  lancT=arg[3]
  lancQ=arg[4]
  lancR=arg[5]
  lancP=arg[6]
  lancB=arg[7]
  lanct=arg[8]
  lancq=arg[9]
  lancr=arg[10]
  lancp=arg[11]
  bla_to_print=arg[12]
  len_lancpos=arg[13]
  title=arg[14]
  if len(arg)> 15:
    extra= arg[15]
  
  par_max_expression_per_gene=0.2 # see parameters.c
  #tr_Rr=[];fr_Rr=[]
  #tr_Rp=[];fr_Rp=[]
  #tr_T=[];fr_T=[]
  #tr_Q=[];fr_Q=[]
  
  # 							Expression_Stuff( lk0,                   lA,               lk1,               lhealth,   par_max_expression_per_gene, lancAct ,lancInact  )
  fr_Rr, lfactinact_Rr, tr_Rr = Expression_Stuff( bla_to_print[4] , bla_to_print[13] , bla_to_print[5], bla_to_print[11], par_max_expression_per_gene, lancR , lancr  )
  fr_Rp, lfactinact_Rp, tr_Rp = Expression_Stuff( bla_to_print[6] , bla_to_print[13] , bla_to_print[7], bla_to_print[11], par_max_expression_per_gene, lancP , lancp  )
  fr_T, lfactinact_T, tr_T = Expression_Stuff( bla_to_print[0] , bla_to_print[13] , bla_to_print[1], bla_to_print[11], par_max_expression_per_gene, lancT , lanct  )
  fr_Q, lfactinact_Q, tr_Q = Expression_Stuff( bla_to_print[2] , bla_to_print[13] , bla_to_print[3], bla_to_print[11], par_max_expression_per_gene, lancQ , lancq  )
  
  #sys.exit(1)
  
  
  #nsubplots= 6 if len(arg)<= 14 else 7
  nsubplots= 7
  fig,ax = plt.subplots(nsubplots, sharex=True)
  
  subplot=0
  #if len(lhist1)>0:
  #  masked_hist = np.ma.masked_where(lhist1 == 0, lhist1)
  #  cmap=truncate_colormap(plt.cm.YlOrRd, minval=0.3, maxval=1., n=100)
  #  print masked_hist.shape
    #cmap=plt.cm.YlOrRd
  #  imgen=ax[subplot].imshow(masked_hist, cmap=cmap, extent=[min(ltime) - 0.5*(ltime[1]-ltime[0])  , max(ltime) + 0.5*(ltime[1]-ltime[0]), 0, nbins], aspect='auto', interpolation='none',label='bla')
  #else:
  #  ax[subplot].plot([1],[1])
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lanc, label='gen size anc',color=color1)
    
  subplot=1
  aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[12], label='S',color=color1)
  aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[13], label='A',color=color2)
  aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[15], label='life span',color=color3)
  aplot=ax[subplot].fill_between( ltime[:len_lancpos],0, bla_to_print[15], label='life span',color=color3,alpha=0.2)
  
  subplot=2
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lancT, label='number of T',color=colorT)
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lancQ, label='number of Q',color=colorQ)
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lancR, label='number of R',color=colorRr)
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lancP, label='number of P',color=colorRp)
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lancB, label='number of B')
  
  subplot=3
  aplot=ax[subplot].plot( ltime[:len_lancpos], lanct, label='number of t',color=colorT)
  aplot=ax[subplot].plot( ltime[:len_lancpos], lancq, label='number of q',color=colorQ)
  aplot=ax[subplot].plot( ltime[:len_lancpos], lancr, label='number of r',color=colorRr)
  aplot=ax[subplot].plot( ltime[:len_lancpos], lancp, label='number of p',color=colorRp)
  
  subplot=4
  aplot=ax[subplot].plot( ltime[:len_lancpos], fr_T, label='fr_T',color=colorT)
  aplot=ax[subplot].plot( ltime[:len_lancpos], fr_Q, label='fr_Q',color=colorQ)
  aplot=ax[subplot].plot( ltime[:len_lancpos], fr_Rr, label='fr_Rr',color=colorRr)
  aplot=ax[subplot].plot( ltime[:len_lancpos], fr_Rp, label='fr_Rp',color=colorRp)
  
  aplot=ax[subplot].plot( ltime[:len_lancpos], lfactinact_T, label='fr act T', linestyle='dotted',color=colorT)
  aplot=ax[subplot].plot( ltime[:len_lancpos], lfactinact_Q, label='fr act Q', linestyle='dotted',color=colorQ)
  aplot=ax[subplot].plot( ltime[:len_lancpos], lfactinact_Rr, label='fr act Rr', linestyle='dotted',color=colorRr)
  aplot=ax[subplot].plot( ltime[:len_lancpos], lfactinact_Rp, label='fr act Rp', linestyle='dotted',color=colorRp)
  
  
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[0], label='kt0')
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[1], label='kt1')
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[2], label='kq0')
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[3], label='kq1')
  
    
  subplot=5
  
  aplot=ax[subplot].plot( ltime[:len_lancpos], tr_T, label='mRNA T',color=colorT)
  aplot=ax[subplot].plot( ltime[:len_lancpos], tr_Q, label='mRNA Q',color=colorQ)
  aplot=ax[subplot].plot( ltime[:len_lancpos], tr_Rr, label='mRNA Rr',color=colorRr)
  aplot=ax[subplot].plot( ltime[:len_lancpos], tr_Rp, label='mRNA Rp',color=colorRp)
  
    
  #aplot=ax[subplot].plot( ltime[:len_lancpos], [ x[0]+x[1]*x[2] for x in zip(bla_to_print[4],bla_to_print[13],bla_to_print[5]) ] , label='Req Rr')
  #aplot=ax[subplot].plot( ltime[:len_lancpos], [  par_max_expression_per_gene * (x[0]+x[1]) for x in zip( lancR,lancr )  ] , label='Max Rr')
  
  #aplot=ax[subplot].plot( ltime[:len_lancpos], fr_Rr, label='tr Rr')
  
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[4], label='kr0')
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[5], label='kr1')
  
  #aplot=ax[subplot].plot( ltime[:len_lancpos], [ x[0]+x[1]*x[2] for x in zip(bla_to_print[6],bla_to_print[13],bla_to_print[7]) ] , label='Req Rp')
  #aplot=ax[subplot].plot( ltime[:len_lancpos], [  par_max_expression_per_gene * (x[0]+x[1]) for x in zip( lancP,lancp )  ] , label='Max Rr')
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[6], label='kp0')
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[7], label='kp1')
  
  subplot=6
  aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[8], label='dupl')
  aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[9], label='del')
  aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[10], label='inact')
  
  aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[11], label='health')
  aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[14], label='av trm')
  
  #subplot with environment taken from parMicro.txt
  #if nsubplots==7:
  #  subplot=6
  #  extra_maxtime=extra[0].index(float(ltime[len_lancpos-1]))
  #  xltime=extra[0][:extra_maxtime]
  #  xenv=extra[1][:extra_maxtime]
  #  aplot=ax[subplot].plot( xltime,xenv, label='S -(!)')
    
  for n in range(nsubplots):
    # Shrink current axis by 20%
    box = ax[n].get_position()
    ax[n].set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax[n].legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize="x-small",ncol=2) # Put a legend to the right of the current axis
    
    cur_xlim = ax[n].get_xlim()	#to set the x axis to start from 0 
    ax[n].set_xlim([min(ltime),cur_xlim[1]])
    
  plt.suptitle(title,fontsize=18)
  
  print len(ltime)
  #plt.xticks(plt.xticks()[0], [ "$"+str(int(t/1000))+"\cdotp 10^3$" if i%2==0 else '' for i,t in enumerate(plt.xticks()[0]) ],size=24)
  #plt.xticks( masked_hist , [ str(time) for time in ltime ])
  #plt.xticks(ltime)
  #plt.legend(loc='upper left')
  plt.savefig("histogenome.png")
  plt.savefig("histogenome.svg")
  plt.show()


def JustPrintData():
  #we use as ttle the name of the executable for the simulation
  files = [f for f in os.listdir('./')]
  title=[f for f in files if 'pMicro' in f][0]
  lin=[]
  with open(anc_data_file, 'r') as fin:
    for line in fin.readlines():
      line=line.split(' ')
      current_time=int(line[0])
      if current_time<minTime: continue
      if current_time>maxTime: break
      lin.append( line )
  
  toplot=zip(*lin)
  ltime=[ int(x) for x in toplot[0]]
  lgenome=toplot[1]
  lanc=[len(genome) for genome in lgenome]
  for genome in lgenome:
    lancB.append( genome.count('B') )
    lancT.append( genome.count('T') )
    lancQ.append( genome.count('Q') )
    lancR.append( genome.count('R') )
    lancP.append( genome.count('P') )
    
    lanct.append( genome.count('t') )
    lancq.append( genome.count('q') )
    lancr.append( genome.count('r') )
    lancp.append( genome.count('p') )
  
  bla_to_print=toplot[2:]
  bla_to_print=[ [ float(bla) for bla in blas ] for blas in bla_to_print ]
  lhist1=[]
  
  #lin2=[]
  
  # NOT needed anymore because we save S (and A) in anc_data_file
  #with open("environment.txt", 'r') as fin:
  #  for line in fin.readlines():
  #    line=line.split(' ')
  #    if float(line[0])>maxTime: break
  #    lin2.append( line )
  
  #if len(lin2)>0:
  #  extra_to_plot=zip(*lin2)
  #  extra_to_plot=[ [float(x) for x in each] for each in extra_to_plot ]
  #else:
  #  print
  #  print "Warning, environment.txt not found, proceeding without"
  #  print
  
  PlotStuff(ltime,lhist1, lanc,lancT,lancQ,lancR,lancP,lancB,lanct,lancq,lancr,lancp,bla_to_print,len(ltime),title)
  sys.exit(0)

    ##############################
#---##########   MAIN   ##########---#
    ##############################
minTime=0 # ONLY USED FOR PRINTING, NOT FOR ANC TRACE
maxTime=999999999
d1={}
ltime=[]
MAXGENOMESIZE=5000

backupdir='backupMicro'
#go through files in current dir
anc_data_file="ancestors_genomes.txt"

if len(sys.argv)>1:
  for i,thisarg in enumerate(sys.argv):
    if thisarg=='-maxtime':
      maxTime=int(sys.argv[i+1])
    if thisarg=='-mintime':
      minTime=int(sys.argv[i+1])
  JustPrintData() #if all you want is printing ancestors_genome.txt

# If the file already exist, see where it stopped last time, 
#  and go on in the backup files from THAT time point
# If the file does not exist, then we start from zero
if os.path.isfile(anc_data_file): 
  print "Output file exists, it is going to be updated"
  bla = subprocess.Popen(['tail', '-n', '1', anc_data_file], stdout=subprocess.PIPE)
  bla = bla.communicate()[0]
  if len(bla)>0:
    initT=int(float(bla.split(' ')[0]))
    update=True
  else:
    initT=0
    update=False
else:
  initT=0
  update=False
##############################
##### --- READ FILES --- #####
##############################
lancdata=[]
files = [f for f in os.listdir('./'+backupdir)]
first_time=True
#backup files are all named backupMicro_t[number]
#we sort this list
for f in sorted(files, key=lambda filename: int( (filename.split('_')[1])[1:] ), reverse=True):
  
  if int( (f.split('_')[1])[1:] ) < initT:
    continue
  
  with open('./'+backupdir+'/'+f,'r') as fp:
    lines=fp.readlines()
    lines=[ line.split(' ') for line in lines ]
    Time = int(lines[0][1])
    
    ltime.append(Time)
    lancdata.append([])
        
    if first_time==True:
      ancpos=[ line[1] for line in lines[1:] ] # we take it everybody, no need to convert to int
      first_time=False
      
    else:
      tags=list(ancpos) #these are the tags that correspond to ancestors in the previous time step (given that we are going backward)
      
      if len(tags)==1:
        # we find the line to which this tag corresponds and saves it, beacause it is ancestral
        
        line=[ x for x in lines if x[0]==tags[0] ][0]
        
        ldata=[ int(line[0]), int(line[1]), line[-1][:-1] ]
        for x in line[2:-1]: ldata.append( float(x) )	# this is eager and takes everything... which can be converted to float
        
        lancdata[-1]=ldata # we save data
        
      #print Time
      #print tags
      ancpos = [ line[1] for line in lines[1:] if line[0] in tags ]
      
#Now we got both ltime and lancdata in reverse order, so we re-reverse them
ltime=ltime[::-1]
lancdata=[x for x in lancdata[::-1] if len(x)!=0]


#print ltime
#print lancdata

#sys.exit(1)
#PRINT the ancestors' genomes as a list 
with open(anc_data_file,'a') as fout:
  if update==True:
    ltime=ltime[1:]	
    lancdata=lancdata[1:]
    #print "Hello",ltime,lancdata
  #else we just print everything
  if len(lancdata)==0:
    print
    print "Nothing to append, this is the end. Goodbye"
    print
    
    sys.exit(0)
  
  for i,ancdata in enumerate(lancdata):
    #print ltime[i], ancdata
    if len(ancdata)==0: break
    genome=ancdata[2]		# position of the genome
    lanc.append(len(genome))
    #print genome, len(genome)
    fout.write(str(ltime[i]))
    for x in ancdata[2:]: 
      fout.write(' '+str(x))
    fout.write('\n')



print 
print "****************************************************************************************"
print 
print "Ancestor tracing complete, re-run the script as ./histo_blabla [-maxtime maxtime] print"
print 
print "****************************************************************************************"
print

