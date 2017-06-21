################################################################################
# Module: test.py
# Description: Test imports and network extraction
# License: GPL3, see full license in LICENSE.txt
# Web: https://github.com/DavidBreuer/CytoSeg
################################################################################

#%%############################################################################# test imports

def test_imports():
    
    import itertools
    import matplotlib.pyplot as plt
    import networkx as nx
    import numpy as np
    import os
    import pandas as pd
    import random
    import scipy as sp
    import scipy.misc
    import scipy.ndimage
    import scipy.optimize
    import scipy.spatial
    import scipy.stats
    import scipy.cluster
    import skimage
    import skimage.filters
    import skimage.io
    import skimage.morphology
    import skimage.feature
    import skimage.segmentation
    import shapely
    import shapely.geometry
    import sys
    import xml
    import xml.dom
    import xml.dom.minidom   
    
    return None

#%%############################################################################# test read tiff
     
def test_read_tiff():
    
    import skimage
    import skimage.io
    im=skimage.io.imread('../examples/data/',plugin='tifffile')  
    
    return None
 
#%%############################################################################# under construction
    
    
   
    



