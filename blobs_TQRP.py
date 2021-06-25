#!/usr/bin/python2.7

'''
Compare different ancestor traces from simulation runs (specified one by one by hand) 
for the last so many time steps (also specified by hand)
'''

import sys,math,os,string
from subprocess import Popen, PIPE
#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np

lanc=[]
lancgen=[]
lancT=[];lancQ=[];lancR=[];lancB=[];lancP=[]
lanct=[];lancq=[];lancr=[];lancp=[]

# COLOR MAP #
#colorT='blue'
#colorQ='green'
#colorRr='red'
#colorRp='cyan'

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


#########################
###                   ###
### --- BLOB PLOT --- ###
###                   ###
#########################
def blob_plot(ax,data,list_pos, colour,datatype,extra=False):
    lextra=[]
    for d,p in zip(data,list_pos):
      
      #this makes the histogram that constitutes the blob plot
      m=min(d) # the largest value, upper edge of the histogram
      M=max(d) # the lowest value, lower...
      if datatype=='int':
        his=np.bincount(d)
        
        his=np.trim_zeros(his)
        his=np.append(his,[0])
        his=np.append([0],his)
        x=np.linspace(m-1,M+1,M-m+3 )
        
        print "his",his
        print "x",x
        '''
        if M-m==0:
          #in these cases histograms are awkard or do not display at all, so we add an extra "zero" data point
          nbins=3
          x = np.linspace(m-1,M+1,nbins) # support for histogram: this is where the histogram is placed
          print "M-m=0; x",x
        elif M-m==1:
          nbins=3
          x = np.linspace(m-1,M,nbins) # support for histogram
          print "M-m=1; x",x
        else:
          nbins=M-m+1
          x = np.linspace(m,M,nbins) # support for histogram
          print "M-m>1; x",x
          
        his,bins = np.histogram(d, bins=x)
        
        #plt.plot(np.append(his,0),bins)
        #plt.show()
      '''
      elif datatype=='float':
        nbins=25
        x = np.linspace(m,M,nbins) # support for histogram of floats
        
        #nbins=25  
        #print "min,max,nbins=",m,M,nbins
      
        his,bins = np.histogram(d, bins=nbins)
      '''  
      if m==M or M-m==1 and datatype=='int': 
        print "m==M or M-m==1"
        print "his",his
        print "bins",bins
      '''     
      
      #scale the histogram to a proper size (you may have to fiddle with this), then place it at position 'pos'
      
      
      #scale_blob=0.5
      
      maxwidht=0.15
      max_his=max(his)
      #print max_his
      if max_his>0.:
        scale_blob=maxwidht/max_his
      else: 
        scale_blob=0.
        print "Warning, max_his is 0." 
        
      shift_his_plus_pos =  [ p + h*scale_blob  for h in his]
      shift_his_minus_pos = [ p - h*scale_blob  for h in his]
      # OLD ->   shift_his_plus_pos =  [ p + h*scale_blob/float(len(d))  for h in his]
      #          shift_his_minus_pos = [ p - h*scale_blob/float(len(d))  for h in his]
      
      color_alpha_blobs=colour
      facecolor,alpha=color_alpha_blobs # assign color and transparency
      #this is the matplotlib function that does the trick
      ax.fill_betweenx(x,shift_his_minus_pos, shift_his_plus_pos, linewidth=0., facecolor= facecolor, alpha=alpha)
      #calculates median or mean, if you want
      if extra=='median': lextra.append( np.median(d) )
      elif extra=='mean': lextra.append( np.mean(d) )
    #and plots it
    if extra != False:
      color = 'orangered'
      ax.plot(list_pos,lextra, color=color_mean_or_median,linestyle='-',marker='D',markersize=5, lw=1.5)
      
 

    ##############################
#---##########   MAIN   ##########---#
    ##############################

d1={}
ltime=[]

howlong=250000
ldir=[]
lhowlong=[]
#go through files in current dir
anc_data_file="ancestors_genomes.txt"
#lifesp_data_file="varMicro.txt"

if len(sys.argv)==1:
  print "Specify some directories and after that for how many generations you want to collect data. default =", howlong 
  print "Also make sure directories exist, and that their path is either in ./ or in ../ or that you are providing a correct absolute path"
elif len(sys.argv)>1:
  #we find path names and howlong
  for x in sys.argv[1:]:
    try:
      bla=int(x)	# if this fails to convert to int then it's probably a dir name
      lhowlong.append(bla)
      howlong=bla
      print "Warning, for now it is not possible to have more than one value for howlong"
    except:
      ldir.append(x)
  
  print lhowlong
  print ldir
  
# Now we check that directories exist
lpath=[]
for d in ldir:
  
  if "/hosts/linuxhome/mutant" in d and os.path.isdir(d): 
    path=d
    if path[-1]!='/': path+="/"	# to make a dir name end with '/'
    lpath.append(path)	# it's an absolute path and we append it as is
    continue
  if "../" in d and os.path.isdir(d):
    print 'hellooooooooooooooooo'
    path=d
    if path[-1]!='/': path+="/"	# to make a dir name end with '/'
    lpath.append(path)	# it's an absolute path and we append it as is
    continue
  
  #alternatively we look for directories either where the program is called 
  # or one level up
  path="./"+d
  if os.path.isdir(path): 
    
    if path[-1]!='/': path+="/"	# to make a dir name end with '/'
    lpath.append(path)
    continue
  
  path="../"+d
  if os.path.isdir(path): 
    
    if path[-1]!='/': path+="/"
    lpath.append(path)
    continue

print "confirmed paths:"
print lpath
print

lfile=[]
lvarMicro_file=[]
# We go get the files
for path in lpath:
  if os.path.isfile(path+anc_data_file): 
    lfile.append(path+anc_data_file)
  else:
    print "ancestor trace file not found in ", path
#  if os.path.isfile(path+lifesp_data_file): 
#    lvarMicro_file.append(path+lifesp_data_file)
#  else:
#   print "varMicro.txt file not found in ", path
  
# and get data (finally!)
lin=[]
ldata=[]
howmanylines=howlong/10 # How long is divided by how frequently ancestors are saved (10 time steps) to give how many lines
						# this is true also for varMicro.txt
minhowmanylines='-'+str(howmanylines)

for bla in lfile:
  lrawdata=[]
  aprocess = Popen(['tail', minhowmanylines, bla], stdout=PIPE, stderr=PIPE)
  stdout, stderr = aprocess.communicate()
  lin=stdout.split('\n')[:-1]
  for line in lin:
    lrawdata.append(line.split(' '))
    #print lrawdata[-1]
  ldata.append(lrawdata)
#print ldata
#sys.exit(1)  
'''
llife=[]
for bla in lvarMicro_file:
  print bla
  lrawdata=[]
  aprocess = Popen(['tail', minhowmanylines, bla], stdout=PIPE, stderr=PIPE)
  stdout, stderr = aprocess.communicate()
  lin=stdout.split('\n')[:-1]
  lrawdata=[ float(line.split(' ')[12]) for line in lin ]
  #print lrawdata[-1]
  llife.append(lrawdata)
'''  
#print llife
#sys.exit(1)  
#PLOTTING

gs = gridspec.GridSpec(3, 3)
ax00 = plt.subplot(gs[0, 0])
ax01 = plt.subplot(gs[0,1])
ax02 = plt.subplot(gs[0,2],sharey=ax01)
ax10 = plt.subplot(gs[1,0])
ax11 = plt.subplot(gs[1,1])
ax12 = plt.subplot(gs[1,2])
axtitle=plt.subplot(gs[2, :])

for i,data in enumerate(ldata):
  #print data
  toplot=zip(*data) 
  #print toplot
  ltime=[ int(x) for x in toplot[0]]
  lgenome=toplot[1]
  lanc=[len(genome) for genome in lgenome]
  
  lancT=[];lancQ=[];lancR=[];lancB=[];lancP=[]
  lanct=[];lancq=[];lancr=[];lancp=[]
  
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
  
  blob_plot(ax00,[lanc],[i],(color1,0.85),'int')
  ax00.set_title('Genome size')
  
  
  print "blobbing T"
  blob_plot(ax01,[lancT],[i-0.3],(colorT,.8),'int')
  print "blobbing Q"
  blob_plot(ax01,[lancQ],[i-0.1],(colorQ,.8),'int')
  print "blobbing Rr"
  blob_plot(ax01,[lancR],[i+0.1],(colorRr,.8),'int')
  print "blobbing Rp"
  blob_plot(ax01,[lancP],[i+0.3],(colorRp,.8),'int')
  ax01.set_title('Active genes')
  
  blob_plot(ax02,[lanct],[i-0.15],(colorT,.8),'int')
  blob_plot(ax02,[lancq],[i-0.05],(colorQ,.8),'int')
  blob_plot(ax02,[lancr],[i+0.05],(colorRr,.8),'int')
  blob_plot(ax02,[lancp],[i+0.15],(colorRp,.8),'int')
  ax02.set_title('Inactive genes')
  
  blob_plot(ax10,[bla_to_print[12]],[i],(color1,.8),'float')
  blob_plot(ax10,[bla_to_print[13]],[i],(color2,.8),'float')
  ax10.set_title('S (blue), A (green) ')
  
  blob_plot(ax11,[bla_to_print[14]],[i],(color1,0.8),'float')
  ax11.set_title('transcr mut.')
  
  blob_plot(ax12, [bla_to_print[15]], [i], (color1,0.8), 'float' )
  ax12.set_title('life span')
  
#ax10.plot(range(10),range(10))
atitle=''
for i,d in enumerate(ldir):
  atitle+=str(i)+': '+d+'\n'
#print suptitle
#plt.suptitle(atitle,fontsize=8);
# place a text box in upper left in axes coords
axtitle.text(0.05, 0.95, '$\\bullet$ T', color= colorT, transform=axtitle.transAxes,fontsize=12, verticalalignment='top')#, bbox=props)
axtitle.text(0.1, 0.95, '$\\bullet$ Q', color= colorQ, transform=axtitle.transAxes,fontsize=12, verticalalignment='top')#, bbox=props)
axtitle.text(0.15, 0.95, '$\\bullet$ Rr', color= colorRr, transform=axtitle.transAxes,fontsize=12, verticalalignment='top')#, bbox=props)
axtitle.text(0.2, 0.95, '$\\bullet$ Rp', color= colorRp, transform=axtitle.transAxes,fontsize=12, verticalalignment='top')#, bbox=props)

axtitle.text(0.05, 0.85, atitle, transform=axtitle.transAxes,fontsize=12, verticalalignment='top')#, bbox=props)

axlim00=ax00.get_ylim()
ax00.set_ylim( 0, axlim00[1])

axlim01=ax01.get_ylim()
ax01.set_ylim( 0, axlim01[1])


axlim11=ax11.get_ylim()
ax11.set_ylim( 0, axlim11[1])
plt.savefig("blobs_TQRP.svg")

plt.show()
#print
#print "NOT SHOWING WHILE DEBUGGING"
#print
print "bye"









































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
  if len(lhist1)>0:
    masked_hist = np.ma.masked_where(lhist1 == 0, lhist1)
    cmap=truncate_colormap(plt.cm.YlOrRd, minval=0.3, maxval=1., n=100)
    print masked_hist.shape
    #cmap=plt.cm.YlOrRd
    imgen=ax[subplot].imshow(masked_hist, cmap=cmap, extent=[min(ltime) - 0.5*(ltime[1]-ltime[0])  , max(ltime) + 0.5*(ltime[1]-ltime[0]), 0, nbins], aspect='auto', interpolation='none',label='bla')
  else:
    ax[subplot].plot([1],[1])
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lanc, label='gen size anc')
  
  subplot=1
  aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[12], label='S')
  aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[13], label='A')
    
  subplot=2
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lancT, label='number of T')
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lancQ, label='number of Q')
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lancR, label='number of R')
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lancP, label='number of P')
  aplot=ax[subplot].plot( ltime[0:len_lancpos], lancB, label='number of B')
  
  subplot=3
  aplot=ax[subplot].plot( ltime[:len_lancpos], lanct, label='number of t')
  aplot=ax[subplot].plot( ltime[:len_lancpos], lancq, label='number of q')
  aplot=ax[subplot].plot( ltime[:len_lancpos], lancr, label='number of r')
  aplot=ax[subplot].plot( ltime[:len_lancpos], lancp, label='number of p')
  
  subplot=4
  aplot=ax[subplot].plot( ltime[:len_lancpos], fr_T, label='fr_T')
  aplot=ax[subplot].plot( ltime[:len_lancpos], fr_Q, label='fr_Q')
  aplot=ax[subplot].plot( ltime[:len_lancpos], fr_Rr, label='fr_Rr')
  aplot=ax[subplot].plot( ltime[:len_lancpos], fr_Rp, label='fr_Rp')
  
  aplot=ax[subplot].plot( ltime[:len_lancpos], lfactinact_T, label='fr act T', linestyle='dotted',color=colorT)
  aplot=ax[subplot].plot( ltime[:len_lancpos], lfactinact_Q, label='fr act Q', linestyle='dotted',color=colorQ)
  aplot=ax[subplot].plot( ltime[:len_lancpos], lfactinact_Rr, label='fr act Rr', linestyle='dotted',color=colorRr)
  aplot=ax[subplot].plot( ltime[:len_lancpos], lfactinact_Rp, label='fr act Rp', linestyle='dotted',color=colorRp)
  
  
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[0], label='kt0')
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[1], label='kt1')
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[2], label='kq0')
  #aplot=ax[subplot].plot( ltime[:len_lancpos], bla_to_print[3], label='kq1')
  
    
  subplot=5
  
  aplot=ax[subplot].plot( ltime[:len_lancpos], tr_T, label='mRNA T')
  aplot=ax[subplot].plot( ltime[:len_lancpos], tr_Q, label='mRNA Q')
  aplot=ax[subplot].plot( ltime[:len_lancpos], tr_Rr, label='mRNA Rr')
  aplot=ax[subplot].plot( ltime[:len_lancpos], tr_Rp, label='mRNA Rp')
  
    
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
    ax[n].set_xlim([0,cur_xlim[1]])
    
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
  title=[f for f in files if 'pMicro0.' in f][0]
  lin=[]
  with open(anc_data_file, 'r') as fin:
    for line in fin.readlines():
      line=line.split(' ')
      if int(line[0])>maxTime: break
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









