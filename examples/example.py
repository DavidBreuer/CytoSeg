################################################################################
# Module: example.py
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
small=25.0                                                                      # smallest components size threshold
factr=0.5                                                                       # fraction of average intensity threshold
depth=7.75                                                                      # spacing between z-slices in xy-pixels spacings (1mum / 0.129mum/pixel = 7.75 pixels)

path='data/'                                                                    # directory with actin and Golgi images                                                        
aa='actin_filter.tif'                                                           # name of actin image 
gg='golgi_filter.tif'                                                           # name of Golgi image

#%%############################################################################# extract and randomize networks

imO=skimage.io.imread(path+aa,plugin='tifffile')[0]                             # open first frame of actin image
imT=skimage.io.imread(path+gg,plugin='tifffile')[0]                             # open Golgi image
mask=skimage.io.imread(path+'mask.tif',plugin='tifffile')>0                     # open mask

track=utils.xmlread(path+'track.xml')                                           # read Golgi tracking results
T=len(track)                                                                    # get number of tracks

print('extract','segment')  
imI=utils.im2d3d(imO)                                                           # if 2D image convert to 3D
imG=skimage.filters.gaussian(imI,sigma)                                         # apply Gaussian filter                                    
imR,imA=utils.skeletonize_graph(imG,mask,sigma,block,small,factr)               # filter and skeletonize actin image
imE=utils.node_graph(imA>0,imG)                                                 # detect filaments and network nodes
     
print('extract','graph')     
gBo,pos=utils.make_graph(imE,imG)                                               # construct graph from filament and node image
gBu=utils.unify_graph(gBo)                                                      # project multigraph to simple graph
gBc=utils.connect_graph(gBu,pos,imG)                                            # connect disconnected components of graph
gBx=utils.centralize_graph(gBc)                                                 # compute edge centrality measures
gBn=utils.normalize_graph(gBx)                                                  # normalize total edge capacity to one    
         
#%%############################################################################# plot data

print('export','plot')  

aspect=2.0                                                                      # set aspect ratio                
alpha=1.0                                                                       # set transparency
lw=1.5                                                                          # set line width        

plt.clf()
gs=mpl.gridspec.GridSpec(1,2,width_ratios=[1,1],height_ratios=[1],left=0.01,bottom=0.01,right=0.99,top=0.99,wspace=0.1,hspace=0.1)                                                       
wh=np.array(np.where(mask))[::-1]
axis=np.hstack(zip(np.nanmin(wh,1),np.nanmax(wh,1)))

plt.subplot(gs[0])                                                              # plot actin image and extracted biological network
plt.title('biological\nactin network')
plt.imshow(imO,cmap='Greys',interpolation='nearest',aspect=aspect)              
ec=1.0*np.array([d['capa'] for u,v,d in gB.edges(data=True)])
nx.draw_networkx_edges(gB,pB[:,:2],edge_color=plt.cm.jet(ec/ec.max()),width=lw,alpha=alpha)
plt.axis(axis)
plt.axis('off')

plt.subplot(gs[1])                                                              # plot Golgi image and tracks
plt.title('Golgi tracks')
plt.imshow(imT,cmap='Greys',interpolation='nearest',aspect=aspect)
for ti,t in enumerate(track[::-1]):
    plt.plot(t[:,1],t[:,2],color=plt.cm.jet(1.0*ti/T),lw=lw,alpha=0.5)
plt.axis(axis)
plt.axis('off')
 
plt.savefig(path+'out_plot.pdf')

