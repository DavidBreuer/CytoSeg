################################################################################
# Module: test.py
# Description: Test imports and network extraction
# License: GPL3, see full license in LICENSE.txt
# Web: https://github.com/DavidBreuer/CytoSeg
################################################################################

#%%############################################################################# imports

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import os
import pandas as pd
import random
import scipy as sp
import skimage
import skimage.io
import skimage.filters
import sys

import utils
    
#%%############################################################################# parameters

sigma=2.0                                                                       # tubeness filter width
block=101.0                                                                     # adaptive median filter block size
small=25.0                                                                      # smallest component size threshold
factr=0.5                                                                       # fraction of average intensity threshold

randw=0                                                                         # randomization method (0 = shuffle edge weights only / 1 = shuffle nodes and edges)
randn=20                                                                        # number of randomized networks

depth=7.75                                                                      # spacing between z-slices in xy-pixels spacings (1mum / 0.129mum/pixel = 7.75 pixels)

path=''                                                                         # directory with actin and Golgi images                                                        
aa='actin_filter.tif'                                                           # name of actin image 
gg='golgi_filter.tif'                                                           # name of Golgi image

#%%############################################################################# extract and randomize networks

imO=skimage.io.imread(path+aa,plugin='tifffile')                                # open actin image
I=len(imO)                                                                      # get number of frames

Z=1                                                                             # set number of z-slices    
shape=imO.shape
if(len(shape)>3): 
    Z=shape[3]

imT=skimage.io.imread(path+gg,plugin='tifffile')                                # open Golgi image

track=utils.xmlread(path+'track.xml')                                           # read Golgi tracking results
T=len(track)                                                                    # get number of tracks
     
mask=skimage.io.imread(path+'mask.tif',plugin='tifffile')>0                     # open mask

R=randn                                                                         # rename number of randomized networks

#%%#

dataB=[]                                                                        # empty list for extracted networks and computed properties
dataP=[]
dataR=[]

for i in range(I):                                                          	# for each frame...

    print('extract',i,I,'segment')  
    imI=imO[i]                                                                  # get actin  image
    imI=utils.im2d3d(imI)                                                       # if 2D image convert to 3D
    imG=skimage.filters.gaussian(imI,sigma)                                     # apply Gaussian filter                                    
    imR,imA=utils.skeletonize_graph(imG,mask,sigma,block,small,factr)           # filter and skeletonize actin image
    imE=utils.node_graph(imA>0,imG)                                             # detect filaments and network nodes
             
    print('extract',i,I,'graph')     
    gBo,pos=utils.make_graph(imE,imG)                                           # construct graph from filament and node image
    gBu=utils.unify_graph(gBo)                                                  # project multigraph to simple graph
    gBc=utils.connect_graph(gBu,pos,imG)                                        # connect disconnected components of graph
    gBx=utils.centralize_graph(gBc)                                             # compute edge centrality measures
    gBn=utils.normalize_graph(gBx)                                              # normalize total edge capacity to one    
    quant=utils.compute_graph(gBn,pos,mask)                                     # compute graph properties
    dataB.append([i,gBn,pos,quant])                                             # append data
        
    for r in range(R):        

        print('extract',i,I,'randomize',r,R)  
        gRo,poz=utils.randomize_graph(gBu,pos,mask,planar=1,weights=randw)      # randomize biological network
        gRu=utils.unify_graph(gRo)                                              # project multigraph to simple graph
        gRc=utils.connect_graph(gRu,poz,imG)                                    # connect disconnected components of graph
        gRx=utils.centralize_graph(gRc)                                         # compute edge centrality measures
        gRn=utils.normalize_graph(gRx)                                          # normalize total edge capacity to one        
        quant=utils.compute_graph(gRn,poz,mask)                                 # compute graph properties         
        dataR.append([i,gRn,poz,quant])                                         # append data
    
#%%############################################################################# plot and export data

print('export','plot')  

i=0                                                                             # choose time point for plotting
r=0                                                                             # choose randomized network for plotting
gB,pB=dataB[i*1+0][1],dataB[i*1+0][2]                                           # get data for biological and randomized network
gR,pR=dataR[i*R+r][1],dataR[i*R+r][2]

plt.clf()
gs=mpl.gridspec.GridSpec(1,3,width_ratios=[1,1,1],height_ratios=[1],left=0.01,bottom=0.01,right=0.99,top=0.99,wspace=0.1,hspace=0.1)
aspect=2.0
alpha=1.0
lw=1.5
wh=np.array(np.where(mask))[::-1]
axis=np.hstack(zip(np.nanmin(wh,1),np.nanmax(wh,1)))

plt.subplot(gs[0])                                                              # plot actin image and extracted biological network
plt.title('biological\nactin network')
plt.imshow(imO[i],cmap='Greys',interpolation='nearest',aspect=aspect)              
ec=1.0*np.array([d['capa'] for u,v,d in gB.edges(data=True)])
nx.draw_networkx_edges(gB,pB[:,:2],edge_color=plt.cm.jet(ec/ec.max()),width=lw,alpha=alpha)
plt.axis(axis)
plt.axis('off')

plt.subplot(gs[1])                                                              # plot actin image and randomized network
plt.title('randomized\nactin network')
plt.imshow(imO[i],cmap='Greys',interpolation='nearest',aspect=aspect)         
ec=1.0*np.array([d['capa'] for u,v,d in gR.edges(data=True)])
nx.draw_networkx_edges(gR,pR[:,:2],edge_color=plt.cm.jet(ec/ec.max()),width=lw,alpha=alpha)
plt.axis(axis)
plt.axis('off')

plt.subplot(gs[2])                                                              # plot Golgi image and tracks
plt.title('Golgi tracks')
plt.imshow(imT[i],cmap='Greys',interpolation='nearest',aspect=aspect)
for ti,t in enumerate(track[::-1]):
    plt.plot(t[:,1],t[:,2],color=plt.cm.jet(1.0*ti/T),lw=lw,alpha=0.5)
plt.axis(axis)
plt.axis('off')
 
plt.savefig(path+'out_plot.pdf')

#%%#

print('export','track')  

idt=np.hstack([np.repeat(ti,len(t)) for ti,t in enumerate(track)])              # convert Golgi tracks to list
tracka=np.vstack([idt,np.vstack(track).T]).T
columns=['ID','t0','x0','y0','z0','avg.intensity0','tot.intensity0','quality0','diameter0','t1','x1','y1','z1','avg.intensity1','tot.intensity1','quality1','diameter1'] # name of recorded Golgi features

df=pd.DataFrame(tracka,columns=columns)                                         # save Golgi tracks as list
df.to_csv(path+'out_track.csv',sep=';',encoding='utf-8')

#%%#

print('export','graph')  

nx.write_gml(gBn,path+'out_graph.gml')                                          # save examplary actin network as graph
 
#%%#

print('export','data')  

quants=['time','# nodes','# edges','# connected components','avg. edge capacity','assortativity','avg. path length','CV path length','algebraic connectivity','CV edge angles','crossing number'] # list of computed network properties

quanta=np.array([np.hstack([d[0],d[-1]]) for d in dataB])                       # save properties of biological networks
df=pd.DataFrame(quanta,columns=quants)
df.to_csv(path+'out_biol.csv',sep=';',encoding='utf-8') 
 
quanta=np.array([np.hstack([d[0],d[-1]]) for d in dataR])                       # save properties of randomized networks
df=pd.DataFrame(quanta,columns=quants)
df.to_csv(path+'out_rand.csv',sep=';',encoding='utf-8') 


