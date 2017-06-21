################################################################################
# Module: fiji.py
# Description: Test imports and network extraction
# License: GPL3, see full license in LICENSE.txt
# Web: https://github.com/DavidBreuer/CytoSeg
################################################################################

#%%############################################################################# imports

import os
import sys
import time

import ij
import ij.IJ
import ij.process
import ij.macro
import ij.measure

import java.io.File

import fiji.plugin.trackmate.Settings as Settings
import fiji.plugin.trackmate.Model as Model
import fiji.plugin.trackmate.TrackMate as TrackMate
import fiji.plugin.trackmate.detection.DogDetectorFactory as DogDetectorFactory
import fiji.plugin.trackmate.tracking.sparselap.SparseLAPTrackerFactory as SparseLAPTrackerFactory
import fiji.plugin.trackmate.tracking.LAPUtils as LAPUtils
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
import fiji.plugin.trackmate.io.TmXmlWriter as TmXmlWriter
import fiji.plugin.trackmate.features.spot.SpotRadiusEstimatorFactory as SpotRadiusEstimatorFactory
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzerFactory as SpotIntensityAnalyzerFactory
import fiji.plugin.trackmate.features.spot.SpotContrastAndSNRAnalyzerFactory as SpotContrastAndSNRAnalyzerFactory
import fiji.plugin.trackmate.features.edges.EdgeTargetAnalyzer as EdgeTargetAnalyzer
import fiji.plugin.trackmate.features.edges.EdgeVelocityAnalyzer as EdgeVelocityAnalyzer
import fiji.plugin.trackmate.features.track.TrackLocationAnalyzer as TrackLocationAnalyzer
import fiji.plugin.trackmate.features.track.TrackIndexAnalyzer as TrackIndexAnalyzer
import fiji.plugin.trackmate.features.track.TrackDurationAnalyzer as TrackDurationAnalyzer
import fiji.plugin.trackmate.features.track.TrackSpeedStatisticsAnalyzer as TrackSpeedStatisticsAnalyzer

#%%############################################################################# parameters

ij.IJ.run('Colors...', 'foreground=white background=black selection=magenta')   # set colors

overwrite_register=1                                                            # overwrite existing files (0 = no / 1 = yes)
overwrite_mask=1
overwrite_filter=1
overwrite_track=1

roll=50                                                                         # rolling ball radius in pixels

radius=4.0                                                                      # Golgi radius in pixels                                                            
quality=4.0                                                                     # minimum quality to accept detects blobs as Golgi
depth=7.75                                                                      # spacing between z-slices in xy-pixels spacings (1mum / 0.129mum/pixel = 7.75 pixels)

distL=24.0                                                                      # maximum linkage distance in pixels
distF=24.0                                                                      # maximum gap-closing distance in pixels
distG=5                                                                         # maximum frame gap number in pixels

path=ij.IJ.getDirectory("Choose a folder")                                      # select directory with actin and Golgi images  
aa='actin_original.tif'                                                         # name of actin image (must contain "original")
gg='golgi_original.tif'                                                         # name of Golgi image (must contain "original")

ij.IJ.run('Close All') 															# close open windows

#%%############################################################################# original -> register

if(os.path.isfile(path+gg.replace('original','register')) and overwrite_register==0): # skip step if output already exists
    print 'register','skip'
    
else:
    print 'register','load'
    for li,l in enumerate([aa,gg]):                                             # for actin and Golgi images...
        ij.IJ.open(path+l)                                                      # open image    
        ij.IJ.run('8-bit')                                                      # convert to 8-bit image
        #ij.IJ.run('Subtract Background...', 'rolling='+str(roll)+' stack')
        #ij.IJ.run('Bleach Correction', 'correction=[Simple Ratio] background=0')
    
    print 'register','merge'   
    ij.IJ.run('Merge Channels...', 'c1='+aa+' c2='+gg)                  		# merge actin and Golgi images
    
    print 'register','split'
    ij.IJ.selectWindow('RGB')													# select merged image
    ij.IJ.run('StackReg', 'transformation=[Rigid Body]')                        # register merged imgage
    ij.IJ.run('Split Channels')                                                 # split image channels
    
    print 'register','save'
    for li,l in enumerate([aa,gg]):                                             # for actin and Golgi images...
        if(li==0): ij.IJ.selectWindow('RGB (red)')
        if(li==1): ij.IJ.selectWindow('RGB (green)')
        ij.IJ.run('8-bit')                                                      # convert to 8-bit image
        ij.IJ.saveAs('Tiff',path+l.replace('original','register'))              # save registered image
    
    print 'register','close'
    #ij.IJ.selectWindow('Log')
    #ij.IJ.run('Close')
    ij.IJ.run('Close All')                                                      # close open windows

#%%############################################################################# register -> mask

if(os.path.isfile(path+'mask.tif') and overwrite_mask==0):                         # skip step if output already exists
    print 'mask','skip'
    
else:
    print 'mask','load'
    for li,l in enumerate([aa,gg]):                                             # for actin and Golgi images...
        ij.IJ.open(path+l.replace('original','register'))                       # open image  
        imp=ij.IJ.getImage()                                                    # rename image  
        imp.setTitle('img'+str(li))
        lx,ly,lc,lz,lt=imp.getDimensions()                                      # get image dimensions
        lz,lt=min(lz,lt),max(lz,lt)
        I=lc*lz*lt
        ij.IJ.run('8-bit')                                                      # convert to 8-bit image
        ij.IJ.run('Properties...', 'channels=1 slices=1 frames='+str(I)+' unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000') # concatenate z-slices and channels as time points (frames)
    
    print 'mask','concatenate'
    ij.IJ.run('Concatenate...', '  title=[Concatenated Stacks] image1=img0 image2=img1 image3=[-- None --]') # concatenate actin and Golgi stacks
    ij.IJ.run('Z Project...', 'projection=[Max Intensity]')                     # z-project combined stack
    
    print 'mask','select'
    ij.IJ.setTool('Freehand')  
    ij.IJ.runMacro("waitForUser ('Pause', 'Please select object of interest.')") # manually select cellular region of interest    
    
    print 'mask','clear'
    ij.IJ.run('Clear Outside')                                                  # clear outside selection
    ij.IJ.run('Fill', 'slice')                                                  # set selected region to uniform intensity
    ij.IJ.run('Measure')                                                        # measure intensity of selected region
    rt=ij.measure.ResultsTable.getResultsTable() 
    ma=rt.getValue('Mean', 0)
    ij.IJ.run('Divide...', 'value='+str(ma))                                    # normalize intensity of selection region to one
    
    print 'mask','save' 
    ij.IJ.run('8-bit')                                                          # convert to 8-bit image      
    ij.IJ.saveAs('Tiff',path+'mask.tif')                                        # save mask
    
    print 'mask','close'
    ij.IJ.selectWindow('Results')                                               # close open windows
    ij.IJ.run('Close')
    ij.IJ.run('Close All')

#%%############################################################################# register -> filter

if(os.path.isfile(path+gg.replace('original','filter')) and overwrite_filter==0): # skip step if output already exists
    print 'filter','skip'
    
else:    
    for li,l in enumerate([aa,gg]):                                             # for actin and Golgi images...  
        print 'filter','load'
        ij.IJ.open(path+l.replace('original','register'))                       # open image  
        imp=ij.IJ.getImage()                                                    # rename image
        imp.setTitle('image')
        lx,ly,lc,lz,lt=imp.getDimensions()                                      # get image dimensions
        lz,lt=min(lz,lt),max(lz,lt)
        I=lc*lz*lt

        print 'filter','filter'
        ij.IJ.run('8-bit')                                                      # convert to 8-bit image
        ij.IJ.run('Properties...', 'channels=1 slices=1 frames='+str(I)+' unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000') # concatenate z-slices and channels as time points (frames) 
        ij.IJ.run('Bleach Correction', 'correction=[Simple Ratio] background=0') # correct photobleaching        
        ij.IJ.run('Subtract Background...', 'rolling='+str(roll)+' stack')      # perform rolling ball background subtraction
        #ij.IJ.run('Despeckle', 'stack')                                        # apply median filter
        #ij.IJ.run('Enhance Contrast...', 'saturated=0.1 process_all use')      # enhance contrast

        print 'filter','crop'
        ij.IJ.open(path+'mask.tif')                                             # open mask
        imp=ij.IJ.getImage()                                                    # rename mask
        imp.setTitle('mask')
        ij.IJ.run('Calculator Plus', 'i1=image i2=mask operation=[Multiply: i2 = (i1*i2) x k1 + k2] k1=1 k2=0 create') # multiply image and mask
        ij.IJ.run('Stack to Hyperstack...', 'order=xyczt(default) channels=1 slices='+str(lz)+' frames='+str(lt)+' display=Grayscale') # convert to hyperstack with lz z-slices and lt time points (frames)

        print 'filter','close'  
        ij.IJ.saveAs('Tiff',path+l.replace('original','filter'))                # save filtered images
        ij.IJ.selectWindow('Log')                                               # close open windows
        ij.IJ.run('Close')
        ij.IJ.run('Close All')

#%%############################################################################# filter -> track     

if(os.path.isfile(path+'track.xml') and overwrite_track==0):                    # skip step if output already exists
    print 'track','skip'
    
else:
    print 'track','load'
    imp=ij.IJ.openImage(path+gg.replace('original','filter'))                   # open Golgi image  
    imp.show()
    lx,ly,lc,lz,lt=imp.getDimensions()                                          # get image dimensions
    lz,lt=min(lz,lt),max(lz,lt)
    ij.IJ.run('Properties...', 'channels=1 slices='+str(lz)+' frames='+str(lt)+' unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth='+str(depth)) # set number of z-slices and frames and spacing between z-slices

    print 'track','setup'
    model=Model()                                                               # set tracking model factories
    settings=Settings()
    settings.setFrom(imp)
    settings.detectorFactory=DogDetectorFactory() 
    settings.detectorSettings['DO_SUBPIXEL_LOCALIZATION']=True                  # do subpixel localization
    settings.detectorSettings['RADIUS']=radius                                  # set blob radius
    settings.detectorSettings['TARGET_CHANNEL']=1                               # set target channel
    settings.detectorSettings['THRESHOLD']=1.0                                  # set detection threshold to one to exclude region outside mask
    settings.detectorSettings['DO_MEDIAN_FILTERING']=True                       # do median filtering
    settings.addSpotAnalyzerFactory(SpotIntensityAnalyzerFactory())
    settings.addSpotAnalyzerFactory(SpotRadiusEstimatorFactory())
    settings.addSpotAnalyzerFactory(SpotContrastAndSNRAnalyzerFactory())    
    settings.initialSpotFilterValue=quality
    settings.addSpotFilter(FeatureFilter('QUALITY',quality,True))  
    settings.trackerFactory=SparseLAPTrackerFactory()
    settings.trackerSettings=LAPUtils.getDefaultLAPSettingsMap()    
    settings.addEdgeAnalyzer(EdgeTargetAnalyzer())
    settings.addEdgeAnalyzer(EdgeVelocityAnalyzer())
    settings.addTrackAnalyzer(TrackSpeedStatisticsAnalyzer())
    settings.addTrackAnalyzer(TrackDurationAnalyzer())
    settings.addTrackAnalyzer(TrackIndexAnalyzer())
    settings.addTrackAnalyzer(TrackLocationAnalyzer())    
    settings.trackerSettings['LINKING_MAX_DISTANCE']=distL                      # set maximum linkage 
    settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE']=distF                  # maximum gap-closing distance in pixels
    settings.trackerSettings['MAX_FRAME_GAP']=distG                             # maximum frame gap number in pixels

    print 'track','tracking'
    trackmate=TrackMate(model,settings)                                         # detect and track Golgi
    trackmate.process()   

    print 'track','save'
    out=java.io.File(path+'track.xml')                                          # save tracking results
    xw=TmXmlWriter(out)
    xw.appendModel(model)
    xw.appendSettings(settings)
    xw.writeToFile()
    ij.IJ.run('Close All')                                                      # close open windows


