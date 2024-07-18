# # ran 1
# import math as math
# import numpy as np
# #import dot2 as dot
# #from py_mentat import *
# #from py_post import *
# import os
# import subprocess
# import pandas as pd
# import csv
# import sys
# import time
# import remoting
# remoting.launchApplication()  
# import apex
# from apex.construct import Point3D,  Point2D
# from pathlib import Path
# import shutil

# #from sympy import *

# #from functions import material2d, material3d
# #from File_paths2 import filepaths

# SELECT THE PARAMETER TO BE CREATED FROM LATIN HYPERCUBE
# The blocks for only two parameter latin hypercube is inside:
# C:\Users\18507522\FinalSBSpring\CODE\3DModelSave\3BlockActuatorFinal
Sel_lat =99

filetree = pd.read_csv("D:\SB_SIMULATIONS\FinalResults\Chamb3_TwoLat\latindata.txt", sep=',',header=None)
Len_A = filetree[0][Sel_lat-1]
Wid_A = filetree[1][Sel_lat-1]
print(Len_A)
print(Wid_A)

# remoting.shutdownApplication()

# deletes default model created incase
model_1 = apex.currentModel()
model_1.close()

# ran 2
# note this seciton of code creates a folder called CODE to contain all subsequent code (if not already created)
import os

os.getcwd()
ob_path = os.getcwd()
# This section of code is to 
# -------------create CODE -------------------
if not os.path.isdir(ob_path+"\\CODE"):
    
    # if the demo_folder2 directory is 
    # not present then create it.
    os.makedirs(ob_path+"\\CODE")
    
root_dir = ob_path+"\\CODE"

#--------------create CODE/FinalOpt ----------------
if not os.path.isdir(root_dir+"\\FinalOpt"):
    
    # if the demo_folder2 directory is 
    # not present then create it.
    os.makedirs(root_dir+"\\FinalOpt")
    
#--------------create CODE/FinalOpt/Results_Opt  ----------------
if not os.path.isdir(root_dir+"\\FinalOpt\\Results_Opt"):
    
    # if the demo_folder2 directory is 
    # not present then create it.
    os.makedirs(root_dir+"\\FinalOpt\\Results_Opt")
    
#--------------part creates the MARC folder for proc folder --------------
if not os.path.isdir(root_dir+"\\FinalOpt\\Marc"):
    
    # if the demo_folder2 directory is 
    # not present then create it.
    os.makedirs(root_dir+"\\FinalOpt\\Marc")
	
	
# ran 3
# NOTE THIS WAS TO TRY AND DELETE:
# creates folder if it doesnt exist then deletes content inside
if not os.path.isdir(root_dir+"\\APEX"):
    
    # if the demo_folder2 directory is 
    # not present then create it.
    os.makedirs(root_dir+"\\APEX")
    
folder = root_dir + "\\APEX"
for filename in os.listdir(folder):
    file_path = os.path.join(folder, filename)
    try:
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)
    except Exception as e:
        print('Failed to delete %s. Reason: %s' % (file_path, e))

print(root_dir)

# # Highlihgting
# from notebook.services.config import ConfigManager
# cm = ConfigManager()
# cm.update('notebook', {'highlight_selected_word': {
#     'delay': 500,
#     'code_cells_only': True,
# }})

selg = 6

f_names = ["_DPLP","_DMLP","_Lm","_DPLM","_Dm","","_DMM","_DMLM","_LMM"]
file = pd.read_csv(root_dir+"\\"+"Fing_points"+f_names[selg-1]+".txt", sep=',',header=None)
file2 = pd.read_csv(root_dir+"\\"+"Bone_points"+f_names[selg-1]+".txt", sep=',',header=None)
LHC = pd.read_csv(root_dir+"\\"+"Latin_Hypercube_Dot.txt", sep=',',header=None, engine='python')

SL0 = file[1][0]

SL1 = file[1][3]
SL2 = file[1][6]
SL3 = file[1][9]
SL4 = file[1][15]

TL0 = -SL0
TL1 = -SL1
TL2 = -SL2
TL3 = -SL3
TL4 = -SL4

# ["_DPLP","_DMLP","_Lm","_DPLM","_Dm","","_DMM","_DMLM","_LMM"]
L_matlab_1 = [0.02168,0.02168,0.01594,0.01693,0.0193,0.0193,0.0193,0.01693,0.02267]
L_matlab_2 = [0.08769,0.08769,0.06694,0.07042,0.07929,0.07929,0.07929,0.07042,0.0907]
L_matlab_3 = [0.04853,0.04853,0.03537,0.03763,0.04308,0.04308,0.04308,0.03763,0.0508]
L_matlab_4 = [0.03612,0.03612,0.02596,0.02769,0.03191,0.03191,0.03191,0.02769,0.03785]

# want to extend the fingers a little bit
L_ext = (-SL4)-(-SL3)

#L0 = round(TL1 - (TL1-TL0)/0.5054   ,7)#0.7454#TL4 +L_matlab_1[selg-1]*1.1725*10.4/25 - L_arc_tot
#L1 = round(TL0 +(TL1-TL0)/2    ,7)
#L2 = round(TL1 +(TL2-TL1)/2    ,7)
#L3 = round(TL2 +(TL3-TL2)/2    ,7)
#L4 = round(TL4 + L_matlab_1[selg-1]*1.1725    ,7)#L_arc4

file1 = open(root_dir+"\\"+"gap_height.txt","r")
gap_height = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"meshing_size.txt","r")
meshing_size = float(file1.read())
file1.close()

meshing_size = 1e-3

file1 = open(root_dir+"\\"+"ch_var1.txt","r")
ch_var1 = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"ch_var2.txt","r")
ch_var2 = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"ch_var3.txt","r")
ch_var3 = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"paper_width.txt","r")
paper_width = float(file1.read())
file1.close()

apex.setScriptUnitSystem(unitSystemName = r'''m-kg-s-N''')

#file1 = open(root_dir+"\\"+"width.txt","r")
#width = float(file1.read())
#file1.close()

#file1 = open(root_dir+"\\"+"height.txt","r")
#height = float(file1.read())
#file1.close()

#file1 = open(root_dir+"\\"+"bredth.txt","r")
#bredth = float(file1.read())
#file1.close()

file1 = open(root_dir+"\\"+"no_ch.txt","r")
no_ch = int(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"ac_th_up.txt","r")
ac_th_up = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"ac_th_side.txt","r")
ac_th_side = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"ac_th_for.txt","r")
ac_th_for = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"gap.txt","r")
gap = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"air_vent.txt","r")
air_vent = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"gap_height.txt","r")
gap_height = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"ends_ext.txt","r")
ends_ext = float(file1.read())
file1.close()

# file1 = open(root_dir+"\\"+"L0.txt","r")
# SL0 = float(file1.read())
# file1.close()

# file1 = open(root_dir+"\\"+"L1.txt","r")
# SL1 = float(file1.read())
# file1.close()

# file1 = open(root_dir+"\\"+"L2.txt","r")
# SL2 = float(file1.read())
# file1.close()

# file1 = open(root_dir+"\\"+"L3.txt","r")
# SL3 = float(file1.read())
# file1.close()

file1 = open(root_dir+"\\"+"fin_rad.txt","r")
fin_rad = float(file1.read())
file1.close

file1 = open(root_dir+"\\"+"meshsize.txt","r")
meshsize = float(file1.read())
file1.close()

meshsize = 0.001

# DPLPBM  
# DPLPBP 
# DMLPB
samp = 101
#LHC = pd.read_csv(root_dir+"\\"+"LHC_100.txt", sep='	',header=None, engine='python')
LHC = pd.read_csv(root_dir+"\\"+"LHC_100_2.txt", sep='	',header=None, engine='python')

LHC_H = 0.5#LHC[0][samp-1]#LHC[1][samp-1] Depth is height
LHC_B = Len_A#LHC[1][samp-1]#LHC[2][samp-1] L is breadth
LHC_W = Wid_A#LHC[3][samp-1]
height = round( 0.01-0.005+0.01*LHC_H      ,7)
bredth = round(10e-3-5e-3+10e-3*LHC_B,7)##round(0.008-0.003+0.006*LHC_B      ,7)  # this is in direction of length
width = round(0.02-0.01+0.02*LHC_W      ,7) #

print(LHC_H)
print(LHC_B)
print(LHC_W)
 
# breadths Lm  0.00593478      0.00574286     0.005711
# breadth LM  0.00745508 0.008229
print('width:',width)
print('height:',height)
print('bredth:',bredth)

L0 = round(TL1 - (TL1-TL0)/0.5054   ,7)
L1 = round(TL1 - (TL1-TL0)/0.5054+3*bredth+gap*3 ,7)

#width: 0.02
#height: 0.01
#bredth: 0.008

# reading additional variables for setting rigid links:
file1 = open(root_dir+"\\"+"rigid1.txt","r")
rigid1 = float(file1.read())
file1.close()
file1 = open(root_dir+"\\"+"rigid2.txt","r")
rigid2 = float(file1.read())
file1.close()
file1 = open(root_dir+"\\"+"rigid3.txt","r")
rigid3 = float(file1.read())
file1.close()

def getAllParts():
    # get and return all solids from current model
    model_1 = apex.currentModel()
    parts = model_1.getParts(True)
    return parts

def getAllSolids():
    solids = apex.entityCollection()
    parts = getAllParts()
    for part in parts:
        solids += part.getSolids()
    return solids

def getAllSurfs():
    surfs = apex.entityCollection()
    parts = getAllParts()
    for part in parts:
        surfs += part.getSurfaces()
    return surfs

# added
def getAllSurfFaces():
    surffaces = apex.entityCollection()
    surfs = getAllSurfs()
    for surf in surfs:
        surffaces += surf.getFaces()
    return surffaces

#
def getAllFaces():
    faces = apex.entityCollection()
    surfs = getAllSurfs()
    for surf in surfs:
        faces += surf.getFaces()
    solids = getAllSolids()
    for solid in solids:
        faces += solid.getFaces()
    return faces

allSurfs = getAllSurfs()

for surf in allSurfs:
    print(surf.name)

ch_width_wings = round(width-2*ac_th_side    ,7) # this just gives how wide the actuator is y axis
ch_height = round(height-ac_th_up    ,7)
ch_bredth = round(bredth-2*ac_th_for    ,7)
#del ch_width

# ==============================================
# This is forcing the number of chambers to be 5
# ==============================================

ch_var1 = 3
ch_var2 = 7
ch_var3 = 5
ch_var4 = 8

lr = L1-L0
chn = ch_var1
tr = lr/chn - 2*ac_th_for-gap
ch_var1 = chn-1+(ac_th_for+tr)/(2*ac_th_for+tr+gap)
tr4 = tr

# lr = L2-L1
# chn = ch_var2
# tr = lr/chn - 2*ac_th_for-gap
# ch_var2 = chn-1+(ac_th_for+tr)/(2*ac_th_for+tr+gap)
# tr4 = tr

# lr = L3-L2
# chn = ch_var3
# tr = lr/chn - 2*ac_th_for-gap
# ch_var3 = chn-1+(ac_th_for+tr)/(2*ac_th_for+tr+gap)
# tr4 = tr

# lr = L4-L3
# chn = ch_var4
# tr = lr/chn - 2*ac_th_for-gap
# ch_var4 = chn-1+(ac_th_for+tr)/(2*ac_th_for+tr+gap)
# tr4 = tr

Origin_start = L0

ch_sel_1 = ch_var1#3.31666666666666665#3.8
no_ch_max_1 = (L1-L0-ac_th_for-gap)/(2*ac_th_for+0.001+gap) # NB CORRECT
if ch_sel_1>no_ch_max_1:
    ch_sel_1 = no_ch_max_1
    print('Dont Proceed')

ch_width_1 = (L1-L0-ac_th_for-gap)/ch_sel_1 - 2*ac_th_for-gap
#print(ch_width_1)

num = math.floor(ch_sel_1)
#print(num)
rem_n = ch_sel_1-num
print(rem_n)
#pointa = ac_th_for/(2*ac_th_for+ch_width_1+gap)
#ointb = (ac_th_for+ch_width_1)/(2*ac_th_for+ch_width_1+gap)
#pointc = (ac_th_for*2+ch_width_1)/(2*ac_th_for+ch_width_1+gap)
#print(pointa)
#print(pointb)
#print(pointc) 

model_1 = apex.currentModel()

# creates blocks without holes part 1
def doSketch_1():
    part_1 = model_1.getCurrentPart()
    if part_1 is None:
        part_1 = model_1.createPart()
    sketch_1 = part_1.createSketchOnGlobalPlane(
        name = 'Sketch 1',
        plane = apex.construct.GlobalPlane.XY,
        alignSketchViewWithViewport = False
    )
    st_pt = round(Origin_start   ,7)
    #varb = Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))
    #varc = Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)


    for a in range(math.floor(ch_sel_1)):
        rectangle_1 = sketch_1.createRectangle2Point(
            name = "Rectangle 1",
            location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(a)   ,7), width/2 ),
            diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(a)   ,7), -width/2 )
        )
    if rem_n == 0:
        # problem: ch_sel_1 without -1 goes over the available so 2 ch code will mode to 3 ch
        # sol: trying to split into two cases
        vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + ac_th_for   ,7)
        varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + (2*ac_th_for+ch_width_1)   ,7)

        # for varc we are subtracting from end so no -1
        varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        # for vard going from end - gap and ac_th_side
        vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - ac_th_side-gap   ,7)
        vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7)
        
        rectangle_1 = sketch_1.createRectangle2Point(
            name = "Rectangle 3",
            location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))\
                           +gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7), -width/2 ),
            diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))\
                           +gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)+0.002   ,7), width/2 )
        )
    else:
        # if segment PNA whole number dont want to make A part of ch +1 segment which doesnt exist thus -1
        vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + ac_th_for   ,7)
        varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + (2*ac_th_for+ch_width_1)    ,7)

        # we are doing +1 here
        varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - ac_th_side - gap   ,7)
        vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1))   ,7)
        print('vard:',vard)
        print('varc:',varc)
        print('varb:',varb)
        print('L1-0.003:',round(L1-gap-ac_th_for,7))
        if rem_n>0:
            if varc <= varb: # very important had to shift this from varc < varb 
                print('varc <= varb')
                rectangle_1 = sketch_1.createRectangle2Point(
                    name = "Rectangle 2",
                    location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7), width/2 ),
                    #diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap) + (2*ac_th_for+ch_width_1+gap)\
                    #                   *(math.floor(ch_sel_1))-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)+0.002   ,7), -width/2 )
                    diagonal = Point2D( round(varc+0.002   ,7), -width/2 )
                )
            if varc > varb:
                print('varc > varb')
                rectangle_1 = sketch_1.createRectangle2Point(
                    name = "Rectangle 2",
                    location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7), width/2 ),
                    diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+\
                                       (2*ac_th_for+ch_width_1)   ,7), -width/2 )
                )

                rectangle_1 = sketch_1.createRectangle2Point(
                    name = "Rectangle 3",
                    location = Point2D( round(varc   ,7), -width/2 ),
                    diagonal = Point2D( round(varc+0.002   ,7), width/2 )
                    #location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap) + (2*ac_th_for+ch_width_1+gap)\
                    #                   *(math.floor(ch_sel_1))-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7), -width/2 ),
                    #diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap) + (2*ac_th_for+ch_width_1+gap)\
                    #                   *(math.floor(ch_sel_1))-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)+0.002   ,7), width/2 )
                
                )

    return sketch_1.completeSketch( fillSketches = True )

newbodies = doSketch_1()
#print(type(newbodies))


a = getAllSurfFaces()
#pushpull
pushpull = apex.geometry.pushPull(
    target = a,
    method = apex.geometry.PushPullMethod.Normal,
    behavior = apex.geometry.PushPullBehavior.Extrude,
    removeInnerLoops = False,
    createBooleanUnion = False,
    distance = height,
    direction = [ 0, 0.0, 1.0 ] )
	
	
if rigid1 == 0:
    model_1 = apex.currentModel()


    def doSketch_1():
        part_1 = model_1.getCurrentPart()
        if part_1 is None:
            part_1 = model_1.createPart()
        sketch_1 = part_1.createSketchOnGlobalPlane(
            name = 'Sketch 1',
            plane = apex.construct.GlobalPlane.XY,
            alignSketchViewWithViewport = False
        )
        st_pt = round(Origin_start  ,7)
        for a in range(math.floor(ch_sel_1)):
            rectangle_1 = sketch_1.createRectangle2Point(
                name = "Rectangle 1",
                location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(a)+ac_th_for   ,7), width/2-ac_th_side ),
                diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(a)-ac_th_for   ,7), -width/2 +ac_th_side)
            )
        if rem_n ==0:
            vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + ac_th_for   ,7)
            varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + (2*ac_th_for+ch_width_1)   ,7)

            # for varc we are subtracting from end so no -1
            varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
            # for vard going from end - gap and ac_th_side
            vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - ac_th_side-gap   ,7)
            vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7)
        else:
            vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + ac_th_for   ,7)
            varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + (2*ac_th_for+ch_width_1)    ,7)

            # we are doing +1 here
            
            varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
            vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - ac_th_side - gap   ,7)
            vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1))   ,7)
            
            print('vard:',vard)
            print('varc:',varc)
            print('varb:',varb)
            print('L1-0.003:', round(L1-gap-ac_th_for,7))
            #=========================REMOVE==================

            #===========================ENDREMOVE=========-==
            #varc > vara
            #varc > vard
            #varc <= varb

            if varc <= vara:
                print('varc <= vara')
            elif varc > vara:
                print('varc > vara')
                if varc <= vard:
                    print('varc <= vard')
                    
                    rectangle_1 = sketch_1.createRectangle2Point(
                        name = "Rectangle 2",
                        location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for   ,7), width/2-ac_th_side ),
                        diagonal = Point2D( round(varc   ,7), -width/2+ac_th_side )
                    )
                elif varc > vard:
                    print('varc > vard')
                    if varc <= varb:
                        print('varc <= varb')
                        rectangle_1 = sketch_1.createRectangle2Point(
                            name = "Rectangle 2",
                            location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for   ,7), width/2-ac_th_side ),
                            #diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)\
                            #                   *(math.floor(ch_sel_1))-ac_th_for   ,7), -width/2+ac_th_side )
                            diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for+ch_width_1\
                                                  ,7), -width/2+ac_th_side )
                        )
                    elif varc > varb:
                        print('varc > varb')
                        rectangle_1 = sketch_1.createRectangle2Point(
                            name = "Rectangle 2",
                            location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for   ,7), width/2-ac_th_side ),
                            diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for+ch_width_1\
                                                  ,7), -width/2+ac_th_side )
                        )

    #     if rem_n == 0:
    #     vara = st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for

    #     varb = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))
    #     varc = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)
    #     vard = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_side
    #     vare = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap
    #     if varc <= vara:
    #         print('nothing')
    #     elif varc > vara:
    #         if varc <= vard:
    #             rectangle_1 = sketch_1.createRectangle2Point(
    #                 name = "Rectangle 2",
    #                 location = Point2D( st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for, width/2-ac_th_side ),
    #                 diagonal = Point2D( varc, -width/2+ac_th_side )
    #             )
    #         elif varc > vard:
    #             if varc <= varb:
    #                 rectangle_1 = sketch_1.createRectangle2Point(
    #                     name = "Rectangle 2",
    #                     location = Point2D( st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for, width/2-ac_th_side ),
    #                     diagonal = Point2D( st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for, -width/2+ac_th_side )
    #                 )
    #             elif varc > varb:
    #                 rectangle_1 = sketch_1.createRectangle2Point(
    #                     name = "Rectangle 2",
    #                     location = Point2D( st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for, width/2-ac_th_side ),
    #                     diagonal = Point2D( st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for, -width/2+ac_th_side )
    #                 )

        return sketch_1.completeSketch( fillSketches = True )

    newbodies = doSketch_1()
    #print(type(newbodies))


    # In[20]:


    # 6
    # THANK YOU T_____T x2 caves
    a = getAllSurfFaces()
    #pushpull
    pushpull = apex.geometry.pushPull(
        target = a,
        method = apex.geometry.PushPullMethod.Normal,
        behavior = apex.geometry.PushPullBehavior.Extrude,
        removeInnerLoops = False,
        createBooleanUnion = False,
        distance = height-ac_th_up,
        direction = [ 0, 0.0, 1.0 ] )

if rigid1 == 0:
    model_1 = apex.currentModel()

    def doSketch_1():
        part_1 = model_1.getCurrentPart()
        if part_1 is None:
            part_1 = model_1.createPart()
        sketch_1 = part_1.createSketchOnGlobalPlane(
            name = 'Sketch 1',
            plane = apex.construct.GlobalPlane.XY,
            alignSketchViewWithViewport = False
        )
        for a in range(math.floor(ch_sel_1)):
            rectangle_1 = sketch_1.createRectangle2Point(
                name = "Rectangle 1",
                location = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) +(2*ac_th_for+ch_width_1+gap)*(a)   ,7), width/2 ),
                diagonal = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(a)+gap   ,7), -width/2 )
            )
        vara = round(Origin_start+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for   ,7)
        varb = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7)
        varc = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        vard = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_side   ,7)
        vare = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap   ,7)

        if varc > varb:
            rectangle_1 = sketch_1.createRectangle2Point(
                name = "Rectangle 1",
                location = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) +(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7), width/2 ),
                diagonal = Point2D( varc, -width/2 )
            )   



        return sketch_1.completeSketch( fillSketches = True )

    newbodies = doSketch_1()

    a = getAllSurfFaces()
    #pushpull
    pushpull = apex.geometry.pushPull(
        target = a,
        method = apex.geometry.PushPullMethod.Normal,
        behavior = apex.geometry.PushPullBehavior.Extrude,
        removeInnerLoops = False,
        createBooleanUnion = False,
        distance = air_vent+gap_height,
        direction = [ 0, 0.0, 1.0 ] )
    
else:    
    model_1 = apex.currentModel()

    def doSketch_1():
        part_1 = model_1.getCurrentPart()
        if part_1 is None:
            part_1 = model_1.createPart()
        sketch_1 = part_1.createSketchOnGlobalPlane(
            name = 'Sketch 1',
            plane = apex.construct.GlobalPlane.XY,
            alignSketchViewWithViewport = False
        )
        for a in range(math.floor(ch_sel_1)):
            rectangle_1 = sketch_1.createRectangle2Point(
                name = "Rectangle 1",
                location = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) +(2*ac_th_for+ch_width_1+gap)*(a)   ,7), width/2 ),
                diagonal = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(a)+gap   ,7), -width/2 )
            )
        vara = round(Origin_start+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for   ,7)
        varb = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7)
        varc = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        vard = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_side   ,7)
        vare = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap   ,7)

        if varc > varb:
            rectangle_1 = sketch_1.createRectangle2Point(
                name = "Rectangle 1",
                location = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) +(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7), width/2 ),
                diagonal = Point2D( varc, -width/2 )
            )   



        return sketch_1.completeSketch( fillSketches = True )

    newbodies = doSketch_1()

    a = getAllSurfFaces()
    #pushpull
    pushpull = apex.geometry.pushPull(
        target = a,
        method = apex.geometry.PushPullMethod.Normal,
        behavior = apex.geometry.PushPullBehavior.Extrude,
        removeInnerLoops = False,
        createBooleanUnion = False,
        distance = height,
        direction = [ 0, 0.0, 1.0 ] )


# Make this shorter

def doSketch_1():
    part_1 = model_1.getCurrentPart()
    if part_1 is None:
        part_1 = model_1.createPart()
    sketch_1 = part_1.createSketchOnGlobalPlane(
        name = 'Sketch 1',
        plane = apex.construct.GlobalPlane.XY,
        alignSketchViewWithViewport = False
    )
    rectangle_1 = sketch_1.createRectangle2Point(
        name = "Rectangle 1",
        location = Point2D( L0, width/2 ),
        diagonal = Point2D( L1-gap, -width/2 )
    )

    return sketch_1.completeSketch( fillSketches = True )

newbodies = doSketch_1()
print(type(newbodies))


# In[57]:


# 6
# THANK YOU T_____T x2 caves
a = getAllSurfFaces()
#pushpull
pushpull = apex.geometry.pushPull(
    target = a,
    method = apex.geometry.PushPullMethod.Normal,
    behavior = apex.geometry.PushPullBehavior.Extrude,
    removeInnerLoops = False,
    createBooleanUnion = False,
    distance = 2e-3,
    direction = [ 0, 0.0, -1.0 ] )


part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()

part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------

solid_1 = part_2.getSolid( name = sol_name )
solid_1.exportGeometry(
    filename =  root_dir+"\\"+"APEX\\"+"Solid1.x_t",
    cadFormat = apex.geometry.CADFormat.ParasolidText,
    unitSystem = "m",
    stlMergeCells = False,
    exportVirtualFaces = False,
    exportVirtualFaceMethod = apex.VirtualFaceExportMethod.NurbsOnly,
    exportBodyPositionInLocal = False
)

# deletes default model created incase
model_1 = apex.currentModel()
model_1.close()

Origin_start = L0

ch_sel_1 = ch_var1#3.31666666666666665#3.8
no_ch_max_1 = (L1-L0-ac_th_for-gap)/(2*ac_th_for+0.001+gap) # NB CORRECT
if ch_sel_1>no_ch_max_1:
    ch_sel_1 = no_ch_max_1
    print('Dont Proceed')

ch_width_1 = (L1-L0-ac_th_for-gap)/ch_sel_1 - 2*ac_th_for-gap
#print(ch_width_1)

num = math.floor(ch_sel_1)
#print(num)
rem_n = ch_sel_1-num
print(rem_n)
#pointa = ac_th_for/(2*ac_th_for+ch_width_1+gap)
#ointb = (ac_th_for+ch_width_1)/(2*ac_th_for+ch_width_1+gap)
#pointc = (ac_th_for*2+ch_width_1)/(2*ac_th_for+ch_width_1+gap)
#print(pointa)
#print(pointb)
#print(pointc) 

model_1 = apex.currentModel()

# creates blocks without holes part 1
def doSketch_1():
    part_1 = model_1.getCurrentPart()
    if part_1 is None:
        part_1 = model_1.createPart()
    sketch_1 = part_1.createSketchOnGlobalPlane(
        name = 'Sketch 1',
        plane = apex.construct.GlobalPlane.XY,
        alignSketchViewWithViewport = False
    )
    st_pt = round(Origin_start   ,7)
    #varb = Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))
    #varc = Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)


    for a in range(math.floor(ch_sel_1)):
        rectangle_1 = sketch_1.createRectangle2Point(
            name = "Rectangle 1",
            location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(a)   ,7), width/2 ),
            diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(a)   ,7), -width/2 )
        )
    if rem_n == 0:
        # problem: ch_sel_1 without -1 goes over the available so 2 ch code will mode to 3 ch
        # sol: trying to split into two cases
        vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + ac_th_for   ,7)
        varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + (2*ac_th_for+ch_width_1)   ,7)

        # for varc we are subtracting from end so no -1
        varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        # for vard going from end - gap and ac_th_side
        vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - ac_th_side-gap   ,7)
        vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7)
        
        rectangle_1 = sketch_1.createRectangle2Point(
            name = "Rectangle 3",
            location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))\
                           +gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7), -width/2 ),
            diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))\
                           +gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)+0.002   ,7), width/2 )
        )
    else:
        # if segment PNA whole number dont want to make A part of ch +1 segment which doesnt exist thus -1
        vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + ac_th_for   ,7)
        varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + (2*ac_th_for+ch_width_1)    ,7)

        # we are doing +1 here
        varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - ac_th_side - gap   ,7)
        vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1))   ,7)
        print('vard:',vard)
        print('varc:',varc)
        print('varb:',varb)
        print('L1-0.003:',round(L1-gap-ac_th_for,7))
        if rem_n>0:
            if varc <= varb: # very important had to shift this from varc < varb 
                print('varc <= varb')
                rectangle_1 = sketch_1.createRectangle2Point(
                    name = "Rectangle 2",
                    location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7), width/2 ),
                    #diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap) + (2*ac_th_for+ch_width_1+gap)\
                    #                   *(math.floor(ch_sel_1))-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)+0.002   ,7), -width/2 )
                    diagonal = Point2D( round(varc+0.002   ,7), -width/2 )
                )
            if varc > varb:
                print('varc > varb')
                rectangle_1 = sketch_1.createRectangle2Point(
                    name = "Rectangle 2",
                    location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7), width/2 ),
                    diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+\
                                       (2*ac_th_for+ch_width_1)   ,7), -width/2 )
                )

                rectangle_1 = sketch_1.createRectangle2Point(
                    name = "Rectangle 3",
                    location = Point2D( round(varc   ,7), -width/2 ),
                    diagonal = Point2D( round(varc+0.002   ,7), width/2 )
                    #location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap) + (2*ac_th_for+ch_width_1+gap)\
                    #                   *(math.floor(ch_sel_1))-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7), -width/2 ),
                    #diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap) + (2*ac_th_for+ch_width_1+gap)\
                    #                   *(math.floor(ch_sel_1))-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)+0.002   ,7), width/2 )
                
                )

    return sketch_1.completeSketch( fillSketches = True )

newbodies = doSketch_1()
#print(type(newbodies))


a = getAllSurfFaces()
#pushpull
pushpull = apex.geometry.pushPull(
    target = a,
    method = apex.geometry.PushPullMethod.Normal,
    behavior = apex.geometry.PushPullBehavior.Extrude,
    removeInnerLoops = False,
    createBooleanUnion = False,
    distance = -height,
    direction = [ 0, 0.0, 1.0 ] )
	
if rigid1 == 0:
    model_1 = apex.currentModel()


    def doSketch_1():
        part_1 = model_1.getCurrentPart()
        if part_1 is None:
            part_1 = model_1.createPart()
        sketch_1 = part_1.createSketchOnGlobalPlane(
            name = 'Sketch 1',
            plane = apex.construct.GlobalPlane.XY,
            alignSketchViewWithViewport = False
        )
        st_pt = round(Origin_start  ,7)
        for a in range(math.floor(ch_sel_1)):
            rectangle_1 = sketch_1.createRectangle2Point(
                name = "Rectangle 1",
                location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(a)+ac_th_for   ,7), width/2-ac_th_side ),
                diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(a)-ac_th_for   ,7), -width/2 +ac_th_side)
            )
        if rem_n ==0:
            vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + ac_th_for   ,7)
            varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + (2*ac_th_for+ch_width_1)   ,7)

            # for varc we are subtracting from end so no -1
            varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
            # for vard going from end - gap and ac_th_side
            vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - ac_th_side-gap   ,7)
            vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7)
        else:
            vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + ac_th_for   ,7)
            varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + (2*ac_th_for+ch_width_1)    ,7)

            # we are doing +1 here
            
            varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
            vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - ac_th_side - gap   ,7)
            vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1))   ,7)
            
            print('vard:',vard)
            print('varc:',varc)
            print('varb:',varb)
            print('L1-0.003:', round(L1-gap-ac_th_for,7))
            #=========================REMOVE==================

            #===========================ENDREMOVE=========-==
            #varc > vara
            #varc > vard
            #varc <= varb

            if varc <= vara:
                print('varc <= vara')
            elif varc > vara:
                print('varc > vara')
                if varc <= vard:
                    print('varc <= vard')
                    
                    rectangle_1 = sketch_1.createRectangle2Point(
                        name = "Rectangle 2",
                        location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for   ,7), width/2-ac_th_side ),
                        diagonal = Point2D( round(varc   ,7), -width/2+ac_th_side )
                    )
                elif varc > vard:
                    print('varc > vard')
                    if varc <= varb:
                        print('varc <= varb')
                        rectangle_1 = sketch_1.createRectangle2Point(
                            name = "Rectangle 2",
                            location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for   ,7), width/2-ac_th_side ),
                            #diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)\
                            #                   *(math.floor(ch_sel_1))-ac_th_for   ,7), -width/2+ac_th_side )
                            diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for+ch_width_1\
                                                  ,7), -width/2+ac_th_side )
                        )
                    elif varc > varb:
                        print('varc > varb')
                        rectangle_1 = sketch_1.createRectangle2Point(
                            name = "Rectangle 2",
                            location = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for   ,7), width/2-ac_th_side ),
                            diagonal = Point2D( round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for+ch_width_1\
                                                  ,7), -width/2+ac_th_side )
                        )

    #     if rem_n == 0:
    #     vara = st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for

    #     varb = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))
    #     varc = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)
    #     vard = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_side
    #     vare = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap
    #     if varc <= vara:
    #         print('nothing')
    #     elif varc > vara:
    #         if varc <= vard:
    #             rectangle_1 = sketch_1.createRectangle2Point(
    #                 name = "Rectangle 2",
    #                 location = Point2D( st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for, width/2-ac_th_side ),
    #                 diagonal = Point2D( varc, -width/2+ac_th_side )
    #             )
    #         elif varc > vard:
    #             if varc <= varb:
    #                 rectangle_1 = sketch_1.createRectangle2Point(
    #                     name = "Rectangle 2",
    #                     location = Point2D( st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for, width/2-ac_th_side ),
    #                     diagonal = Point2D( st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for, -width/2+ac_th_side )
    #                 )
    #             elif varc > varb:
    #                 rectangle_1 = sketch_1.createRectangle2Point(
    #                     name = "Rectangle 2",
    #                     location = Point2D( st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for, width/2-ac_th_side ),
    #                     diagonal = Point2D( st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for, -width/2+ac_th_side )
    #                 )

        return sketch_1.completeSketch( fillSketches = True )

    newbodies = doSketch_1()
    #print(type(newbodies))


    # In[20]:


    # 6
    # THANK YOU T_____T x2 caves
    a = getAllSurfFaces()
    #pushpull
    pushpull = apex.geometry.pushPull(
        target = a,
        method = apex.geometry.PushPullMethod.Normal,
        behavior = apex.geometry.PushPullBehavior.Extrude,
        removeInnerLoops = False,
        createBooleanUnion = False,
        distance = -(height-ac_th_up),
        direction = [ 0, 0.0, 1.0 ] )


if rigid1 == 0:
    model_1 = apex.currentModel()

    def doSketch_1():
        part_1 = model_1.getCurrentPart()
        if part_1 is None:
            part_1 = model_1.createPart()
        sketch_1 = part_1.createSketchOnGlobalPlane(
            name = 'Sketch 1',
            plane = apex.construct.GlobalPlane.XY,
            alignSketchViewWithViewport = False
        )
        for a in range(math.floor(ch_sel_1)):
            rectangle_1 = sketch_1.createRectangle2Point(
                name = "Rectangle 1",
                location = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) +(2*ac_th_for+ch_width_1+gap)*(a)   ,7), width/2 ),
                diagonal = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(a)+gap   ,7), -width/2 )
            )
        vara = round(Origin_start+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for   ,7)
        varb = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7)
        varc = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        vard = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_side   ,7)
        vare = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap   ,7)

        if varc > varb:
            rectangle_1 = sketch_1.createRectangle2Point(
                name = "Rectangle 1",
                location = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) +(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7), width/2 ),
                diagonal = Point2D( varc, -width/2 )
            )   



        return sketch_1.completeSketch( fillSketches = True )

    newbodies = doSketch_1()

    a = getAllSurfFaces()
    #pushpull
    pushpull = apex.geometry.pushPull(
        target = a,
        method = apex.geometry.PushPullMethod.Normal,
        behavior = apex.geometry.PushPullBehavior.Extrude,
        removeInnerLoops = False,
        createBooleanUnion = False,
        distance = -(air_vent+gap_height),
        direction = [ 0, 0.0, 1.0 ] )
    
else:    
    model_1 = apex.currentModel()

    def doSketch_1():
        part_1 = model_1.getCurrentPart()
        if part_1 is None:
            part_1 = model_1.createPart()
        sketch_1 = part_1.createSketchOnGlobalPlane(
            name = 'Sketch 1',
            plane = apex.construct.GlobalPlane.XY,
            alignSketchViewWithViewport = False
        )
        for a in range(math.floor(ch_sel_1)):
            rectangle_1 = sketch_1.createRectangle2Point(
                name = "Rectangle 1",
                location = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) +(2*ac_th_for+ch_width_1+gap)*(a)   ,7), width/2 ),
                diagonal = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(a)+gap   ,7), -width/2 )
            )
        vara = round(Origin_start+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for   ,7)
        varb = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7)
        varc = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap-(2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        vard = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_side   ,7)
        vare = round(Origin_start+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+gap   ,7)

        if varc > varb:
            rectangle_1 = sketch_1.createRectangle2Point(
                name = "Rectangle 1",
                location = Point2D( round(Origin_start+(2*ac_th_for+ch_width_1) +(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7), width/2 ),
                diagonal = Point2D( varc, -width/2 )
            )   



        return sketch_1.completeSketch( fillSketches = True )

    newbodies = doSketch_1()

    a = getAllSurfFaces()
    #pushpull
    pushpull = apex.geometry.pushPull(
        target = a,
        method = apex.geometry.PushPullMethod.Normal,
        behavior = apex.geometry.PushPullBehavior.Extrude,
        removeInnerLoops = False,
        createBooleanUnion = False,
        distance = -height,
        direction = [ 0, 0.0, 1.0 ] )



def doSketch_1():
    part_1 = model_1.getCurrentPart()
    if part_1 is None:
        part_1 = model_1.createPart()
    sketch_1 = part_1.createSketchOnGlobalPlane(
        name = 'Sketch 1',
        plane = apex.construct.GlobalPlane.XY,
        alignSketchViewWithViewport = False
    )
    rectangle_1 = sketch_1.createRectangle2Point(
        name = "Rectangle 1",
        location = Point2D( L0, width/2 ),
        diagonal = Point2D( L1-gap, -width/2 )
    )

    return sketch_1.completeSketch( fillSketches = True )

newbodies = doSketch_1()
print(type(newbodies))


# In[57]:


# 6
# THANK YOU T_____T x2 caves
a = getAllSurfFaces()
#pushpull
pushpull = apex.geometry.pushPull(
    target = a,
    method = apex.geometry.PushPullMethod.Normal,
    behavior = apex.geometry.PushPullBehavior.Extrude,
    removeInnerLoops = False,
    createBooleanUnion = False,
    distance = -(2e-3),
    direction = [ 0, 0.0, -1.0 ] )


move_allParts = getAllParts()
move_Part_name1 = str(move_allParts[0].name)

move_allSol = getAllSolids()
move_sol_name1 = str(move_allSol[0].name)



_entities = apex.currentModel().getEntities(pathNames=[move_Part_name1+'/'+move_sol_name1])
newEntities = apex.transformTranslate(
    target = _entities,
    direction = [ 0.0, 0.0, -1.000000000000000 ],
    distance = 4e-3,
    makeCopy = False
)

# _target = apex.EntityCollection()
# model_1 = apex.currentModel()
# part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
# allSol = getAllSolids()
# sol_name = str(allSol[0].name)

# export the paper layer
part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()

part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------

solid_1 = part_2.getSolid( name = sol_name )
solid_1.exportGeometry(
    filename = root_dir+"\\"+"APEX\\"+"Solid2.x_t",
    cadFormat = apex.geometry.CADFormat.ParasolidText,
    unitSystem = "m",
    stlMergeCells = False,
    exportVirtualFaces = False,
    exportVirtualFaceMethod = apex.VirtualFaceExportMethod.NurbsOnly,
    exportBodyPositionInLocal = False
)

# deletes default model created incase
model_1 = apex.currentModel()
model_1.close()

# export the botsil layer
model_1 = apex.currentModel()

_geometryFileNames = [
    root_dir+"\\"+"APEX\\"+"Solid1.x_t",
    root_dir+"\\"+"APEX\\"+"Solid2.x_t"
]
_importFilter = {
    
}
model_1.importGeometry(
    geometryFileNames = _geometryFileNames,
    importSolids = True,     #defaults to True
    importSurfaces = True,     #defaults to True
    importCurves = True,     #defaults to True
    importPoints = True,     #defaults to True
    importGeneralBodies = True,     #defaults to True
    importHiddenGeometry = False,     #defaults to False
    importCoordinate = False,     #defaults to True
    importDatumPlane = False,     #defaults to True
    cleanOnImport = True,     #defaults to True
    removeRedundantTopoOnimport = False,     #defaults to False
    loadCompleteTopology = True,     #defaults to True
    sewOnImport = False,     #defaults to False
    skipUnmodified = False,     #defaults to False
    importFilter = _importFilter,
    preview = False,     #defaults to False
    importReviewMode3dxmlCleanOnImport = True,     #defaults to True
    importReviewMode3dxmlSplitOnFeatureVertexAngle = True,     #defaults to True
    importReviewMode3dxmlFeatureAngle = 4.000000000000000e+01,     #defaults to 40.0 degrees
    importReviewMode3dxmlVertexAngle = 4.000000000000000e+01,     #defaults to 40.0 degrees
    importReviewMode3dxmlDetectMachinedFaces = False,     #defaults to False
    importAttributes = False,     #defaults to False
    importPublications = False     #defaults to False
)

# combines to one part
allSol = getAllSolids()
sol_name1 = str(allSol[0].name)
sol_name2 = str(allSol[1].name)

allParts = getAllParts()
Part_name1 = str(allParts[0].name)
Part_name2 = str(allParts[1].name)

model_1 = apex.currentModel()
solid_1 = apex.getSolid( pathName = model_1.name + "/" + Part_name1 + "/"+ sol_name1)
updatedSolidName_ = solid_1.update( name = "Solid1" )

solid_2 = apex.getSolid( pathName = model_1.name + "/" + Part_name2 + "/"+ sol_name2)
updatedSolidName_ = solid_2.update( name = "Solid2" )

part_1 = apex.getPart( pathName = model_1.name+"/" + Part_name1 )
part_2 = apex.getPart( pathName = model_1.name+"/" + Part_name2 )

# selects solid 2 in part 2 and moves it to first
entities_2 = apex.EntityCollection()
entities_2.append(solid_2)
apex.reparent( parent = part_1, target = entities_2 )

# delete the empty parts
entities_del1 = apex.EntityCollection()
entities_del1.append(part_2)
apex.deleteEntities(entities_del1)

#name the solids again
#-------------------- name
allSol = getAllSolids()
sol_name1 = str(allSol[0].name)
sol_name2 = str(allSol[1].name)

# merge the solids
_target = apex.EntityCollection()
#part_4 = apex.getPart( pathName = apex.currentModel().name + "/Part 1" )
solid_1 = part_1.getSolid( name = sol_name1 )
_target.append( solid_1 )
solid_2 = part_1.getSolid( name = sol_name2 )
_target.append( solid_2 )
result = apex.geometry.mergeBoolean(
    target = _target,
    retainOriginalBodies = False,
    mergeSolidsAsCells = True
)

propertyelement3d_1 = apex.catalog.createPropertiesElement3D(
    name = '3D Element Property 1',
    primaryProperties3D = 'PSOLID',
)

propertyelement3d_1 = apex.catalog.getPropertiesElement3D( name = "3D Element Property 1" )
propertyelement3d_1.update(
    name = 'Actuator',
)

def getAllCells():
    cells = apex.entityCollection()
    solids = getAllSolids()
    for solid in solids:
        cells += solid.getCells()
    return cells

allCells = getAllCells()

# assigning materials:
propertyelement3d_4 = apex.catalog.getPropertiesElement3D( name = "Actuator" )
_target = apex.EntityCollection()
_geoms = apex.EntityCollection()
_geoms.extend( solid_1.getCells( ids = str(allCells[0].id) ) )
_target.extend( _geoms )
apex.attribute.assignPropertiesElement3D(target = _target, property = propertyelement3d_4 )

propertyelement3d_4 = apex.catalog.getPropertiesElement3D( name = "Actuator" )
_target = apex.EntityCollection()
_geoms = apex.EntityCollection()
_geoms.extend( solid_1.getCells( ids = str(allCells[1].id) ) )
_target.extend( _geoms )
apex.attribute.assignPropertiesElement3D(target = _target, property = propertyelement3d_4 )

def doSketch_1():
    part_1 = model_1.getCurrentPart()
    if part_1 is None:
        part_1 = model_1.createPart()
    sketch_1 = part_1.createSketchOnGlobalPlane(
        name = 'Sketch 1',
        plane = apex.construct.GlobalPlane.XY,
        alignSketchViewWithViewport = False
    )
    rectangle_1 = sketch_1.createRectangle2Point(
        name = "Rectangle 1",
        location = Point2D( L0, width/2 ),
        diagonal = Point2D( L1-gap, -width/2 )
    )

    return sketch_1.completeSketch( fillSketches = True )

newbodies = doSketch_1()



move_allSurf = getAllSurfs()
move_surf_name1 = str(move_allSurf[0].name)



_entities = apex.currentModel().getEntities(pathNames=[move_Part_name1+'/'+move_surf_name1])
newEntities = apex.transformTranslate(
    target = _entities,
    direction = [ 0.0, 0.0, -1.000000000000000 ],
    distance = 2e-3,
    makeCopy = False
)

# create the paper layer 2D material
propertyelement2d_1 = apex.catalog.createPropertiesElement2D(
    name = '2D Element Property 1',
    primaryProperties2D = 'PSHELL',
)

propertyelement2d_1 = apex.catalog.getPropertiesElement2D( name = "2D Element Property 1" )
propertyelement2d_1.update(
    name = 'Paper',
)

# assign paper layer
allParts = getAllParts()
Part_name1 = str(allParts[0].name)
model_1 = apex.currentModel()

propertyelement2d_2 = apex.catalog.getPropertiesElement2D( name = "Paper" )
_target = apex.EntityCollection()
_geoms = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name + "/"+Part_name1 )
surface_1 = part_1.getSurface( name = move_allSurf[0].name )
_geoms.append( surface_1 )
_target.extend( _geoms )
apex.attribute.assignProperty2D( property2d = propertyelement2d_2, target = _target )

# CREATES GAPS BETWEEN CHAMBERS SPLIT
part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()

part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------
allSurf = getAllSurfs()
surf_name = str(allSurf[0].name)

solid_1 = part_2.getSolid( name = sol_name )
surf_1 = part_2.getSurface( name = surf_name)

# This next part is to get the faces
#--------------------
facecollection = part_1.getSolid( name = sol_name).getFaces()

faceselection = 'Null'

temp_X = 0
temp_Y = 0
temp_Z = 0

id_col = []
x_col = []
y_col = []
z_col = []

iter_count = 0
top_count = 0

dada = []
flag = 0

side1 = width
side2 = height-gap_height-air_vent
for face in facecollection:
    if round(float(face.getArea()),7) == round(side1*side2,7):
        #print(str(round(face.getArea(),4))+"first check " + str(round(height*width,4)))
        #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))
        if round(float(face.getCentroid().getY()),7) == round(0.0,7):
            flag = 0
            for da in range(len(dada)):
                if dada[da] == round(float(face.getCentroid().getX()),7):
                    flag = 1
                    
            if flag != 1:
                dada.append(round(float(face.getCentroid().getX()),7))

                faceselection = str(face.id)
                temp_X = float(face.getCentroid().getX())
                temp_Y = float(face.getCentroid().getY())
                temp_Z = float(face.getCentroid().getZ())
                id_col.append(faceselection)
                x_col.append(temp_X)
                y_col.append(temp_Y)
                z_col.append(temp_Z)
                iter_count += 1
            #print('face selection is: '+ faceselection)

            #--------------------
if faceselection == 'Null':
    raise TypeError("Face wasn't found")

for a in range(iter_count):            
    #print(str(temp_X)+" "+str(temp_Y)+" "+str(temp_Z)+" face: "+str(faceselection))

    _target = apex.EntityCollection()
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    solid_1 = part_1.getSolid( name = sol_name )
    _target.append( solid_1 )
    _location = apex.Coordinate( x_col[a], y_col[a], z_col[a])
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    solid_1 = part_1.getSolid( name = sol_name )
    _snapEntity = solid_1.getFace( id = int(id_col[a]) )
    result = apex.geometry.splitOnFeaturePlane(
        target = _target,
        location = _location,
        snapEntity = _snapEntity,
        splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
    )
    _target = apex.EntityCollection()
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    surf_1 = part_1.getSurface( name = surf_name )
    _target.append( surf_1 )
    result = apex.geometry.splitOnFeaturePlane(
        target = _target,
        location = _location,
        snapEntity = _snapEntity,
        splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
    )

    # gaps split in x direction
part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()
part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------

solid_1 = part_2.getSolid( name = sol_name )

# This next part is to get the faces
#--------------------
facecollection = part_1.getSolid( name = sol_name).getFaces()

faceselection = 'Null'

# splits inside chamber ends
temp_X = 0
temp_Y = 0
temp_Z = 0

id_col = []
x_col = []
y_col = []
z_col = []

iter_count = 0
top_count = 0

dada = []
flag = 0

side1 = ch_width_wings
side2 = height-ac_th_up
for face in facecollection:
    #print('A',round(float(face.getArea()),7))
    #print('B',round(side1*side2,7))
    #print('=================')
    if round(float(face.getArea()),7) == round(side1*side2,7):
        #print('Pass')
        #print(str(round(face.getArea(),4))+"first check " + str(round(height*width,4)))
        #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))
        #print('Centroid',round(float(face.getCentroid().getY()),4))
        if round(float(face.getCentroid().getY()),7) == round(0.0,7):
            
            flag = 0
            for da in range(len(dada)):
                if dada[da] == round(float(face.getCentroid().getX()),7):
                    flag = 1
                    
            if flag != 1:
                dada.append(round(float(face.getCentroid().getX()),7))
                faceselection = str(face.id)
                temp_X = float(face.getCentroid().getX())
                temp_Y = float(face.getCentroid().getY())
                temp_Z = float(face.getCentroid().getZ())
                id_col.append(faceselection)
                x_col.append(temp_X)
                y_col.append(temp_Y)
                z_col.append(temp_Z)
                iter_count += 1
            #print('face selection is: '+ faceselection)

            #--------------------
if faceselection == 'Null':
    raise TypeError("Face wasn't found")

for a in range(iter_count):            
    #print(str(temp_X)+" "+str(temp_Y)+" "+str(temp_Z)+" face: "+str(faceselection))

    _target = apex.EntityCollection()
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    solid_1 = part_1.getSolid( name = sol_name )
    _target.append( solid_1 )
    _location = apex.Coordinate( x_col[a], y_col[a], z_col[a])
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    solid_1 = part_1.getSolid( name = sol_name )
    _snapEntity = solid_1.getFace( id = int(id_col[a]) )
    result = apex.geometry.splitOnFeaturePlane(
        target = _target,
        location = _location,
        snapEntity = _snapEntity,
        splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
    )
    _target = apex.EntityCollection()
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    surf_1 = part_1.getSurface( name = surf_name )
    _target.append( surf_1 )
    result = apex.geometry.splitOnFeaturePlane(
        target = _target,
        location = _location,
        snapEntity = _snapEntity,
        splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
    )
# does nothing?
# Creates splits ends of inside chamber
temp_X = 0
temp_Y = 0
temp_Z = 0

id_col = []
x_col = []
y_col = []
z_col = []

iter_count = 0
top_count = 0

side1 = ch_width_wings
#side2 = height-air_vent
side2 = ch_height
for face in facecollection:
    if round(float(face.getArea()),7) == round(side1*side2-air_vent*air_vent,7):
        #print(str(round(face.getArea(),4))+"first check " + str(round(height*width,4)))
        #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))
        if round(float(face.getCentroid().getY()),4) == round(0.0,4):
            print('obtained:'+str(round(float(face.getArea()),7)))
            print('desired:'+str(round(side1*side2,7)))
            print('')
            #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))    
            faceselection = str(face.id)
            temp_X = float(face.getCentroid().getX())
            temp_Y = float(face.getCentroid().getY())
            temp_Z = float(face.getCentroid().getZ())
            id_col.append(faceselection)
            x_col.append(temp_X)
            y_col.append(temp_Y)
            z_col.append(temp_Z)
            iter_count += 1
            #print('face selection is: '+ faceselection)

            #--------------------
if faceselection == 'Null':
    raise TypeError("Face wasn't found")

for a in range(iter_count):            
    print(str(temp_X)+" "+str(temp_Y)+" "+str(temp_Z)+" face: "+str(faceselection))

    _target = apex.EntityCollection()
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    solid_1 = part_1.getSolid( name = sol_name )
    _target.append( solid_1 )
    _location = apex.Coordinate( x_col[a], y_col[a], z_col[a])
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    solid_1 = part_1.getSolid( name = sol_name )
    _snapEntity = solid_1.getFace( id = int(id_col[a]) )
    result = apex.geometry.splitOnFeaturePlane(
        target = _target,
        location = _location,
        snapEntity = _snapEntity,
        splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
    )

# gaps split in x direction
part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()
part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------

solid_1 = part_2.getSolid( name = sol_name )

# This next part is to get the faces
#--------------------
facecollection = part_1.getSolid( name = sol_name).getFaces()

faceselection = 'Null'

#******************************top split ********************x direct
# This code is used to create the splits needed for meshing
part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()
part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------

solid_1 = part_2.getSolid( name = sol_name )

# This next part is to get the faces
#--------------------
facecollection = part_1.getSolid( name = sol_name).getFaces()

faceselection = 'Null'
temp_X = 0
temp_Y = 0
temp_Z = 0

ch_sel_1 = ch_var1#3.31666666666666665#3.8
no_ch_max_1 = (L1-L0-ac_th_for-gap)/(2*ac_th_for+0.001+gap) # NB CORRECT
if ch_sel_1>no_ch_max_1:
    ch_sel_1 = no_ch_max_1
    print('Dont Proceed')

ch_width_1 = (L1-L0-ac_th_for-gap)/ch_sel_1 - 2*ac_th_for-gap

count = 0
for face in facecollection:
    #print(count)
    #print(str(round((width-2*ac_th_side)*ch_width_1,7)))
    #print("found:"+str(round(float(face.getArea()),7) ))
    if round(float(face.getArea()),7) == round((width-2*ac_th_side)*ch_width_1,7):
        #print(str(round(face.getArea(),8))+"first check " + str(round(height*width,8)))
        #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))
        if round(count,5) == round(0,5):
            if round(float(face.getCentroid().getZ()),7) == round(height-ac_th_up,7):
                count = 1
                #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))    
                faceselection = str(face.id)
                temp_X = float(face.getCentroid().getX())
                temp_Y = float(face.getCentroid().getY())
                temp_Z = float(face.getCentroid().getZ())
                #print('face selection is: '+ faceselection)

            #--------------------
if faceselection == 'Null':
    raise TypeError("Face wasn't found")

print(str(temp_X)+" "+str(temp_Y)+" "+str(temp_Z)+" face: "+str(faceselection))

_target = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
solid_1 = part_1.getSolid( name = sol_name )
_target.append( solid_1 )
_location = apex.Coordinate( temp_X, temp_Y, temp_Z )
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
solid_1 = part_1.getSolid( name = sol_name )
_snapEntity = solid_1.getFace( id = int(faceselection) )
result = apex.geometry.splitOnFeaturePlane(
    target = _target,
    location = _location,
    snapEntity = _snapEntity,
    splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
)
# Doing this again as it has a bottom reflection side as well
#******************************top split ********************x direct
# This code is used to create the splits needed for meshing

flag_g = 0

part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()
part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------

solid_1 = part_2.getSolid( name = sol_name )

# This next part is to get the faces
#--------------------
facecollection = part_1.getSolid( name = sol_name).getFaces()

faceselection = 'Null'
temp_X = 0
temp_Y = 0
temp_Z = 0

ch_sel_1 = ch_var1#3.31666666666666665#3.8
no_ch_max_1 = (L1-L0-ac_th_for-gap)/(2*ac_th_for+0.001+gap) # NB CORRECT
if ch_sel_1>no_ch_max_1:
    ch_sel_1 = no_ch_max_1
    print('Dont Proceed')

ch_width_1 = (L1-L0-ac_th_for-gap)/ch_sel_1 - 2*ac_th_for-gap

count = 0
for face in facecollection:
    #print(count)
    #print(str(round((width-2*ac_th_side)*ch_width_1,7)))
    #print("found:"+str(round(float(face.getArea()),7) ))
    if round(float(face.getArea()),7) == round((width-2*ac_th_side)*ch_width_1,7):
        #print(str(round(face.getArea(),8))+"first check " + str(round(height*width,8)))
        #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))
        if round(count,5) == round(0,5):
            if round(float(face.getCentroid().getZ()),7) == -0.004-round(height-ac_th_up,7):
                count = 1
                #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))    
                faceselection = str(face.id)
                temp_X = float(face.getCentroid().getX())
                temp_Y = float(face.getCentroid().getY())
                temp_Z = float(face.getCentroid().getZ())
                #print('face selection is: '+ faceselection)

            #--------------------
if faceselection == 'Null':
    #raise TypeError("Face wasn't found")
    flag_g = 1

if flag_g == 0:
    print(str(temp_X)+" "+str(temp_Y)+" "+str(temp_Z)+" face: "+str(faceselection))

    _target = apex.EntityCollection()
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    solid_1 = part_1.getSolid( name = sol_name )
    _target.append( solid_1 )
    _location = apex.Coordinate( temp_X, temp_Y, temp_Z )
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    solid_1 = part_1.getSolid( name = sol_name )
    _snapEntity = solid_1.getFace( id = int(faceselection) )
    result = apex.geometry.splitOnFeaturePlane(
        target = _target,
        location = _location,
        snapEntity = _snapEntity,
        splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
    )

#******************************top split ********************at gap
# This code is used to create the splits needed for meshing
part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()
part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------

solid_1 = part_2.getSolid( name = sol_name )

# This next part is to get the faces
#--------------------
facecollection = part_1.getSolid( name = sol_name).getFaces()

faceselection = 'Null'
temp_X = 0
temp_Y = 0
temp_Z = 0

ch_sel_1 = ch_var1#3.31666666666666665#3.8
no_ch_max_1 = (L1-L0-ac_th_for-gap)/(2*ac_th_for+0.001+gap) # NB CORRECT
if ch_sel_1>no_ch_max_1:
    ch_sel_1 = no_ch_max_1
    print('Dont Proceed')

ch_width_1 = (L1-L0-ac_th_for-gap)/ch_sel_1 - 2*ac_th_for-gap

count = 0
for face in facecollection:
    #print(count)
    #print(str(round((width-2*ac_th_side)*ch_width_1,7)))
    #print("found:"+str(round(float(face.getArea()),7) ))
    if round(float(face.getArea()),7) == round((width)*gap,7):
        #print(str(round(face.getArea(),8))+"first check " + str(round(height*width,8)))
        #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))
        if round(count,5) == round(0,5):
            if round(float(face.getCentroid().getZ()),7) == round(air_vent+gap_height,7):
                count = 1
                #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))    
                faceselection = str(face.id)
                temp_X = float(face.getCentroid().getX())
                temp_Y = float(face.getCentroid().getY())
                temp_Z = float(face.getCentroid().getZ())
                #print('face selection is: '+ faceselection)

            #--------------------
if faceselection == 'Null':
    raise TypeError("Face wasn't found")

print(str(temp_X)+" "+str(temp_Y)+" "+str(temp_Z)+" face: "+str(faceselection))

_target = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
solid_1 = part_1.getSolid( name = sol_name )
_target.append( solid_1 )
_location = apex.Coordinate( temp_X, temp_Y, temp_Z )
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
solid_1 = part_1.getSolid( name = sol_name )
_snapEntity = solid_1.getFace( id = int(faceselection) )
result = apex.geometry.splitOnFeaturePlane(
    target = _target,
    location = _location,
    snapEntity = _snapEntity,
    splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
)
# Doing this again as it has a bottom reflection side as well
#******************************top split ********************x direct
# This code is used to create the splits needed for meshing
part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()
part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------

solid_1 = part_2.getSolid( name = sol_name )

# This next part is to get the faces
#--------------------
facecollection = part_1.getSolid( name = sol_name).getFaces()

faceselection = 'Null'
temp_X = 0
temp_Y = 0
temp_Z = 0

ch_sel_1 = ch_var1#3.31666666666666665#3.8
no_ch_max_1 = (L1-L0-ac_th_for-gap)/(2*ac_th_for+0.001+gap) # NB CORRECT
if ch_sel_1>no_ch_max_1:
    ch_sel_1 = no_ch_max_1
    print('Dont Proceed')

ch_width_1 = (L1-L0-ac_th_for-gap)/ch_sel_1 - 2*ac_th_for-gap

count = 0
for face in facecollection:
    #print(count)
    #print(str(round((width-2*ac_th_side)*ch_width_1,7)))
    #print("found:"+str(round(float(face.getArea()),7) ))
    if round(float(face.getArea()),7) == round((width)*gap,7):
        #print(str(round(face.getArea(),8))+"first check " + str(round(height*width,8)))
        #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))
        if round(count,5) == round(0,5):
            if round(float(face.getCentroid().getZ()),7) == -0.004-round(air_vent+gap_height,7):
                count = 1
                #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))    
                faceselection = str(face.id)
                temp_X = float(face.getCentroid().getX())
                temp_Y = float(face.getCentroid().getY())
                temp_Z = float(face.getCentroid().getZ())
                #print('face selection is: '+ faceselection)

            #--------------------
if faceselection == 'Null':
    raise TypeError("Face wasn't found")

print(str(temp_X)+" "+str(temp_Y)+" "+str(temp_Z)+" face: "+str(faceselection))

_target = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
solid_1 = part_1.getSolid( name = sol_name )
_target.append( solid_1 )
_location = apex.Coordinate( temp_X, temp_Y, temp_Z )
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
solid_1 = part_1.getSolid( name = sol_name )
_snapEntity = solid_1.getFace( id = int(faceselection) )
result = apex.geometry.splitOnFeaturePlane(
    target = _target,
    location = _location,
    snapEntity = _snapEntity,
    splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
)

# This code is used to create the splits needed for meshing
part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()
part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------

solid_1 = part_2.getSolid( name = sol_name )

# This next part is to get the faces
#--------------------
facecollection = part_1.getSolid( name = sol_name).getFaces()

faceselection = 'Null'

temp_X = 0
temp_Y = 0
temp_Z = 0

ch_sel_1 = ch_var1#3.31666666666666665#3.8
no_ch_max_1 = (L1-L0-0.002-0.002)/(2*ac_th_for+0.001+gap) # NB CORRECT
if ch_sel_1>no_ch_max_1:
    ch_sel_1 = no_ch_max_1
    
ch_bredth_p = (L1-L0-0.002-0.002)/ch_sel_1 - 2*ac_th_for-gap 

side1 = ch_bredth_p
side2 = height-ac_th_up-air_vent-gap_height
count = 0

for face in facecollection:
    if round(float(face.getArea()),4) == round(side1*side2,4):
        #print(str(round(face.getArea(),4))+"first check " + str(round(height*width,4)))
        #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))
        if count == 0:

            if round(float(face.getCentroid().getY()),4) == round(ch_width_wings/2,4):
                count = 1
                #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))    
                faceselection = str(face.id)
                temp_X = float(face.getCentroid().getX())
                temp_Y = float(face.getCentroid().getY())
                temp_Z = float(face.getCentroid().getZ())
            #print('face selection is: '+ faceselection)

            #--------------------
if faceselection == 'Null':
    raise TypeError("Face wasn't found")

print(str(temp_X)+" "+str(temp_Y)+" "+str(temp_Z)+" face: "+str(faceselection))

_target = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
solid_1 = part_1.getSolid( name = sol_name )
_target.append( solid_1 )
_location = apex.Coordinate( temp_X, temp_Y, temp_Z )
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
solid_1 = part_1.getSolid( name = sol_name )
_snapEntity = solid_1.getFace( id = int(faceselection) )
result = apex.geometry.splitOnFeaturePlane(
    target = _target,
    location = _location,
    snapEntity = _snapEntity,
    splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
)

_target = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
surf_1 = part_1.getSurface( name = surf_name )
_target.append( surf_1 )
result = apex.geometry.splitOnFeaturePlane(
    target = _target,
    location = _location,
    snapEntity = _snapEntity,
    splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
)


# This code is used to create the splits needed for meshing
part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()
part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------

solid_1 = part_2.getSolid( name = sol_name )

# This next part is to get the faces
#--------------------
facecollection = part_1.getSolid( name = sol_name).getFaces()

faceselection = 'Null'

temp_X = 0
temp_Y = 0
temp_Z = 0

ch_sel_1 = ch_var1#3.31666666666666665#3.8
no_ch_max_1 = (L1-L0-0.002-0.002)/(2*ac_th_for+0.001+gap) # NB CORRECT
if ch_sel_1>no_ch_max_1:
    ch_sel_1 = no_ch_max_1
    
ch_bredth_p = (L1-L0-0.002-0.002)/ch_sel_1 - 2*ac_th_for-gap 

side1 = ch_bredth_p
side2 = height-ac_th_up-air_vent-gap_height
count = 0

for face in facecollection:
    if round(float(face.getArea()),4) == round(side1*side2,4):
        #print(str(round(face.getArea(),4))+"first check " + str(round(height*width,4)))
        #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))
        if count == 0:

            if round(float(face.getCentroid().getY()),4) == round(-ch_width_wings/2,4):
                count = 1
                #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))    
                faceselection = str(face.id)
                temp_X = float(face.getCentroid().getX())
                temp_Y = float(face.getCentroid().getY())
                temp_Z = float(face.getCentroid().getZ())
            #print('face selection is: '+ faceselection)

            #--------------------
if faceselection == 'Null':
    raise TypeError("Face wasn't found")

print(str(temp_X)+" "+str(temp_Y)+" "+str(temp_Z)+" face: "+str(faceselection))

_target = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
solid_1 = part_1.getSolid( name = sol_name )
_target.append( solid_1 )
_location = apex.Coordinate( temp_X, temp_Y, temp_Z )
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
solid_1 = part_1.getSolid( name = sol_name )
_snapEntity = solid_1.getFace( id = int(faceselection) )
result = apex.geometry.splitOnFeaturePlane(
    target = _target,
    location = _location,
    snapEntity = _snapEntity,
    splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
)

_target = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
surf_1 = part_1.getSurface( name = surf_name )
_target.append( surf_1 )
result = apex.geometry.splitOnFeaturePlane(
    target = _target,
    location = _location,
    snapEntity = _snapEntity,
    splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
)
#CREATES GAPS BETWEEN CHAMBERS SPLIT
part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()

part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------
allSurf = getAllSurfs()
surf_name = str(allSurf[0].name)

solid_1 = part_2.getSolid( name = sol_name )
surf_1 = part_2.getSurface( name = surf_name)

# This next part is to get the faces
#--------------------
facecollection = part_1.getSolid( name = sol_name).getFaces()

faceselection = 'Null'

temp_X = 0
temp_Y = 0
temp_Z = 0

id_col = []
x_col = []
y_col = []
z_col = []

iter_count = 0
top_count = 0

dada = []
flag = 0

side1 = width-2*ac_th_side
side2 = height-gap_height-air_vent-ac_th_up
for face in facecollection:
    if round(float(face.getArea()),7) == round(side1*side2,7):
        #print(str(round(face.getArea(),4))+"first check " + str(round(height*width,4)))
        #print(str(face.getCentroid().getX())+ " second "+ str(bredth*(a+1)+gap*(a)-bredth/2))
        if round(float(face.getCentroid().getY()),7) == round(0.0,7):
            flag = 0
            for da in range(len(dada)):
                if dada[da] == round(float(face.getCentroid().getX()),7):
                    flag = 1
                    
            if flag != 1:
                if round(face.getCentroid().getX(),7) == round(L0+2*ac_th_for+ch_width_1,7):
                    dada.append(round(float(face.getCentroid().getX()),7))

                    faceselection = str(face.id)
                    temp_X = float(face.getCentroid().getX())
                    temp_Y = float(face.getCentroid().getY())
                    temp_Z = float(face.getCentroid().getZ())
                    id_col.append(faceselection)
                    x_col.append(temp_X)
                    y_col.append(temp_Y)
                    z_col.append(temp_Z)
                    iter_count += 1
                if round(face.getCentroid().getX(),7) == round(L0+2*ac_th_for+ch_width_1+gap+2*ac_th_for+ch_width_1,7):
                    dada.append(round(float(face.getCentroid().getX()),7))

                    faceselection = str(face.id)
                    temp_X = float(face.getCentroid().getX())
                    temp_Y = float(face.getCentroid().getY())
                    temp_Z = float(face.getCentroid().getZ())
                    id_col.append(faceselection)
                    x_col.append(temp_X)
                    y_col.append(temp_Y)
                    z_col.append(temp_Z)
                    iter_count += 1
            #print('face selection is: '+ faceselection)

            #--------------------
if faceselection == 'Null':
    raise TypeError("Face wasn't found")

for a in range(iter_count):            
    #print(str(temp_X)+" "+str(temp_Y)+" "+str(temp_Z)+" face: "+str(faceselection))
    print(str(int(id_col[a])))
    _target = apex.EntityCollection()
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    solid_1 = part_1.getSolid( name = sol_name )
    _target.append( solid_1 )
    _splitter = apex.EntityCollection()
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1")
    solid_1 = part_1.getSolid( name = sol_name )
    _splitter.extend( solid_1.getFaces( ids = str(int(id_col[a])) ) )
    result = apex.geometry.splitEntityWithOffsetFaces(
        target = _target,
        splitter = _splitter,
        offset = -5.000000000000000e-04,
        splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
    )
    print('x:',x_col[a])
    print('ot:',L0+2*ac_th_for+ch_width_1)
	
    _target = apex.EntityCollection()
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
    surf_1 = part_1.getSurface( name = surf_name )
    _target.append( surf_1 )
    _splitter = apex.EntityCollection()
    part_1 = apex.getPart( pathName = model_1.name+"/Part 1")
    solid_1 = part_1.getSolid( name = sol_name )
    _splitter.extend( solid_1.getFaces( ids = str(int(id_col[a])) ) )
    result = apex.geometry.splitEntityWithOffsetFaces(
        target = _target,
        splitter = _splitter,
        offset = -5.000000000000000e-04,
        splitBehavior = apex.geometry.GeometrySplitBehavior.Partition
    )

    # gaps split in x direction
part_1 = model_1.getCurrentPart()
if part_1 is None:
    part_1 = model_1.createPart()
part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

#-------------------- name
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------

solid_1 = part_2.getSolid( name = sol_name )

# This next part is to get the faces
#--------------------
facecollection = part_1.getSolid( name = sol_name).getFaces()

faceselection = 'Null'
#=======================cutting half ====================
#@@@@@@@@@@@@@@@@@@@@@@
x_min =round(L0-2e-3-meshsize/2,7)
x_max = round( L1+2e-3+meshsize/2,7)# act_end
y_min = round(-width/2-meshsize/2,7)
y_max = round(0,6)#str(round(meshsize/2,6))
z_min = round(-height-0.004-meshsize/2,7)
z_max = round(height+meshsize/2,7)

result = apex.geometry.createBoxByLocationOrientation(
    name = 'Box',
    description = '',
    length = x_max-x_min,
    height = y_max-y_min,
    depth = z_max-z_min,
    origin = apex.Coordinate( x_min, y_min, z_min ),
    orientation = apex.Orientation(0.0, 0.0, 0.0)
    )

_target = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
allSol = getAllSolids()
sol_name = str(allSol[0].name)

part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

solid_1 = part_2.getSolid( name = sol_name )
_target.append( solid_1 )
_subtractingEntities = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name + "/Part 1" )
body_1 = part_1.getGeometryBody( name = "Box 1" )
_subtractingEntities.append( body_1 )
result = apex.geometry.subtractBoolean(
    target = _target,
    subtractingEntity = _subtractingEntities,
    retainOriginalBodies = True
)


_target = apex.EntityCollection()
part_1 = apex.getPart( pathName = apex.currentModel().name + "/Part 1" )
surface_1 = part_1.getSurface( name = surf_name )
_target.append( surf_1 )
_subtractingEntities = apex.EntityCollection()
part_1 = apex.getPart( pathName = apex.currentModel().name + "/Part 1" )
body_1 = part_1.getGeometryBody( name = "Box 1" )
_subtractingEntities.append( body_1 )
result = apex.geometry.subtractBoolean(
    target = _target,
    subtractingEntity = _subtractingEntities,
    retainOriginalBodies = True
)


geomes_1 = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name + "/Part 1" )
body_1 = part_1.getGeometryBody( name = "Box 1" )
geomes_1.append( body_1 )
apex.deleteEntities(geomes_1)

#@@@@@@@@@@@@@@@@@@@@@@
#========================================================







#meshing_size = 0.001

#=======================================================
#*****************************************************

#                 CREATE MESH AND EXPORT

#***********************************************************
#========================================================

_target = apex.EntityCollection()
part_1 = apex.getPart( pathName = model_1.name+"/Part 1" )
allSol = getAllSolids()
sol_name = str(allSol[0].name)
#--------------------
part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

solid_1 = part_2.getSolid( name = sol_name )

_target.append( solid_1 )
_SweepFace = apex.EntityCollection()
result = apex.mesh.createHexMesh(
    name = "",
    target = _target,
    meshSize = meshing_size,
    surfaceMeshMethod = apex.mesh.SurfaceMeshMethod.Auto,
    mappedMeshDominanceLevel = 2,
    elementOrder = apex.mesh.ElementOrder.Linear,
    refineMeshUsingCurvature = False,
    elementGeometryDeviationRatio = 0.10,
    elementMinEdgeLengthRatio = 0.20,
    createFeatureMeshOnWashers = False,
    createFeatureMeshOnArbitraryHoles = False,
    preserveWasherThroughMesh = True,
    sweepFace = _SweepFace,
    hexMeshMethod = apex.mesh.HexMeshMethod.Auto,
    projectMidsideNodesToGeometry = True
)

_target = apex.EntityCollection()
allSurf = getAllSurfs()
surf_name = str(allSurf[0].name)
print(surf_name)
part_2 = apex.getPart( pathName = model_1.name+"/Part 1" )

surface_1 = part_2.getSurface( name = surf_name )
_target.append( surface_1 )
featuremeshtypes_1 = apex.mesh.FeatureMeshTypeVector()
result = apex.mesh.createSurfaceMesh(
    name = "",
    target = _target,
    meshSize = meshing_size,
    meshType = apex.mesh.SurfaceMeshElementShape.Quadrilateral,
    meshMethod = apex.mesh.SurfaceMeshMethod.Auto,
    mappedMeshDominanceLevel = 2,
    elementOrder = apex.mesh.ElementOrder.Linear,
    allQuadBoundary = False,
    refineMeshUsingCurvature = False,
    curvatureType = apex.mesh.CurvatureType.EdgeOnly,
    elementGeometryDeviationRatio = 0.10,
    elementMinEdgeLengthRatio = 0.20,
    proximityRefinement = False,
    growFaceMeshSize = False,
    faceMeshGrowthRatio = 1.2,
    createFeatureMeshes = False,
    featureMeshTypes = featuremeshtypes_1,
    projectMidsideNodesToGeometry = True,
    useMeshFlowOptimization = False,
    meshFlow = apex.mesh.MeshFlow.Grid,
    minimalMesh = False
)
mesh1 = apex.getMesh( pathName = apex.currentModel().name + "/Part 1/Mesh 1" )
mesh1.exportFEModel(
    filename = "C:/Users/18507522/FinalSBSpring/CODE/3DModelSave/3BlockActuatorFinal/Sel"+str(Sel_lat)+"_3D.bdf",
    unitSystem = "m-kg-s-N-K",
    exportProperty = apex.ExportProperty.AsDefined
 )

mesh1 = apex.getMesh( pathName = apex.currentModel().name + "/Part 1/Mesh 2" )
mesh1.exportFEModel(
    filename = "C:/Users/18507522/FinalSBSpring/CODE/3DModelSave/3BlockActuatorFinal/Sel"+str(Sel_lat)+"_2D.bdf",
    unitSystem = "m-kg-s-N-K",
    exportProperty = apex.ExportProperty.AsDefined
 )
