from py_post import *
from py_mentat import *
import subprocess
import sys
import math
import time
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

### Default Values ####
pressure = float(50000)#float(130000)
max_load_case = 0.5
#######################

os.getcwd()
ob_path = os.getcwd()
root_dir = ob_path+"\\CODE"

# code creates the folder
if not os.path.isdir(ob_path+"\\CODE\\ExampleResultsForMarc"):
    
    # if the demo_folder2 directory is 
    # not present then create it.
    os.makedirs(ob_path+"\\CODE\\ExampleResultsForMarc")
    
root_dir = ob_path+"\\CODE"

file1 = open(root_dir+"\\"+"width.txt","r")
width = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"height.txt","r")
height = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"bredth.txt","r")
bredth = float(file1.read())
file1.close()

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

file1 = open(root_dir+"\\"+"L0.txt","r")
SL0 = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"L1.txt","r")
SL1 = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"L2.txt","r")
SL2 = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"L3.txt","r")
SL3 = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"fin_rad.txt","r")
fin_rad = float(file1.read())
file1.close

file1 = open(root_dir+"\\"+"meshsize.txt","r")
meshsize = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"ch_var1.txt","r")
ch_var1 = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"ch_var2.txt","r")
ch_var2 = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"ch_var3.txt","r")
ch_var3 = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"gap_height.txt","r")
gap_height = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"meshing_size.txt","r")
meshing_size = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"paper_width.txt","r")
paper_width = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"cham_end.txt","r")
cham_end = float(file1.read())
file1.close()

#file1 = open(root_dir+"\\"+"rotstiff.txt","r")
#rotstiff = float(file1.read())
#file1.close()
rotstiff = 0.0001

file1 = open(root_dir+"\\"+"Lm1.txt","r")
Lm1 = float(file1.read())
file1.close()

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

file1 = open(root_dir+"\\"+"quit.txt","r")
quit = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"gen_res.txt","r")
gen_res = float(file1.read())
file1.close()

file1 = open(root_dir+"\\"+"multithread_flag.txt","r")
multithread_flag = float(file1.read())
file1.close()

##### Setting some values ######
# if finger_strap_flag == 1 or 2 with 2 additionally x is also constrained:
strap_offset1 = (SL1-SL0)/2
strap_offset2 = (SL2-SL1)/2
strap_offset3 = (SL3-SL2)/2
FL0= SL0
FL1 = SL1
FL2 = SL2
FL3 = SL3
# Now we are changing the values of L0,L1,L2,L3,Lm1 to consider between joints and not at the joint

Sel_lat =1

filetree = pd.read_csv("D:\SB_SIMULATIONS\FinalResults\Chamb3_TwoLat\latindata.txt", sep=',',header=None)
Len_A = filetree[0][Sel_lat-1]
Wid_A = filetree[1][Sel_lat-1]
print(Len_A)
print(Wid_A)


LHC = pd.read_csv(root_dir+"\\"+"LHC_100_2.txt", sep='	',header=None, engine='python')
LHC_P = 1#LHC[0][samp-1]

LHC_H = 0.5#LHC[0][samp-1]#LHC[1][samp-1] Depth is height
LHC_B = Len_A#1.0#LHC[1][samp-1]#LHC[2][samp-1] L is breadth
LHC_W = Wid_A#0.5#LHC[3][samp-1]

LHC_Var = 0.5
dim = 'H'
if LHC_H != 0.5:
    LHC_Var = LHC_H
    dim = 'H'

if LHC_B != 0.5:
    LHC_Var = LHC_B
    dim = 'B'
    
if LHC_W != 0.5:
    LHC_Var = LHC_W
    dim = 'W'
    
# Dimension to be working on
dim = 'S'
dim_c = Sel_lat#str(int(LHC_Var*100))

#LHC_H = LHC[1][samp-1]
#LHC_B = LHC[2][samp-1]
#LHC_W = LHC[3][samp-1]
height = round( 0.01-0.005+0.01*LHC_H      ,7) #Depth finger
bredth = round(10e-3-5e-3+10e-3*LHC_B,7)# Length               #round(0.008-0.003+0.006*LHC_B      ,7)  # this is in direction of length
width = round(0.02-0.01+0.02*LHC_W      ,7) # Breadth 0.01
#height = round( 0.01-0.005+0.01*LHC_H      ,7) 
#bredth = round(6.875e-3,7)##round(0.008-0.003+0.006*LHC_B      ,7)
#width = round(0.02-0.005+0.01*LHC_W      ,7)


L0 = Lm1
ch_width = width-2*ac_th_side
ch_height = height-ac_th_up
meshsize = 0.001#meshing_size
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

L0 = round(TL1 - (TL1-TL0)/0.5054   ,7)#0.7454#TL4 +L_matlab_1[selg-1]*1.1725*10.4/25 - L_arc_tot
L1 = round(TL0 +(TL1-TL0)/2    ,7)
L2 = round(TL1 +(TL2-TL1)/2    ,7)
L3 = round(TL2 +(TL3-TL2)/2    ,7)
L4 = round(TL4 + L_matlab_1[selg-1]*1.1725    ,7)#L_arc4

Origin_start = L0



#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# DPLPBM  
# DPLPBP 
# DMLPB


 
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
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

ch_width_wings = round(width-2*ac_th_side    ,7) # this just gives how wide the actuator is y axis
ch_height = round(height-ac_th_up    ,7)
ch_bredth = round(bredth-2*ac_th_for    ,7)
#del ch_width

# ==============================================
# This is forcing the number of chambers to be 5
# ==============================================

ch_var1 = 3 # 11
ch_var2 = 7
ch_var3 = 5
ch_var4 = 8

lr = L1-L0
chn = ch_var1
tr = lr/chn - 2*ac_th_for-gap
ch_var1 = chn-1+(ac_th_for+tr)/(2*ac_th_for+tr+gap)
tr4 = tr

lr = L2-L1
chn = ch_var2
tr = lr/chn - 2*ac_th_for-gap
ch_var2 = chn-1+(ac_th_for+tr)/(2*ac_th_for+tr+gap)
tr4 = tr

lr = L3-L2
chn = ch_var3
tr = lr/chn - 2*ac_th_for-gap
ch_var3 = chn-1+(ac_th_for+tr)/(2*ac_th_for+tr+gap)
tr4 = tr

lr = L4-L3
chn = ch_var4
tr = lr/chn - 2*ac_th_for-gap
ch_var4 = chn-1+(ac_th_for+tr)/(2*ac_th_for+tr+gap)
tr4 = tr

Fing_D3_CMC_T = file[2][1]

# THIS IS THE AMOUNT WHICH FINGERS NEED TO BE DISPLACED TO BE UNDER ACTUATOR
move_amount = -Fing_D3_CMC_T-height-2e-3-2e-3

Fing_D3_MCP_C = file[2][3]-Fing_D3_CMC_T-height-2e-3-2e-3
Fing_D3_PIP_C = file[2][6]-Fing_D3_CMC_T-height-2e-3-2e-3
Fing_D3_DIP_C = file[2][9]-Fing_D3_CMC_T-height-2e-3-2e-3

# specify the spring values 
f_names = ["_DPLP","_DMLP","_Lm","_DPLM","_Dm","","_DMM","_DMLM","_LMM"]
r1x = [0.0582,0.085,0.0375,0.043,0.0641,0.0485,0.0465,0.0412,0.0793]
r2x = [0.028,0.0305,0.0196,0.0226,0.0251,0.0252,0.0253,0.0198,0.0309]
r3x = [0.0036,0.0033,0.0023,0.0028,0.0026,0.003,0.0033,0.0021,0.0037]

#1 D+L+,2*D-L+, 3>Lm, 4>D+L-, 5>Dm, 6>O, 7>D-L+, 8>DM, 9>D-L-, 10>LM,  
# SELECT THE VAR USING selg!!!!


mcp_sp = r1x[selg-1]

pip_sp = r2x[selg-1]

dip_sp = r3x[selg-1]

# if set to 1 then will change to including spheres
spheres_flag = 0
# if set to 1 then will consider default change to 1 and change strap_offset for altered
finger_strap_flag = 1

#----------Finger_accurate_flag
#use 1 for with spheres
#use 2 for Hex (2nd order?)
#use 3 for 2nd order Tet
#use 4 for 1st order Tet
finger_accurate_flag = 0 # default 3

# *-* I used 4 - 1

# flag to correct the symmetry in the model
cor_sym_flag = 1

#loadshedding counter
restart_flag = 0

finger_offset = -0.001

# to create the beams set to 1 if set to 0 it will only create the nodes
beams_code_flag = 1

core_thread = 2

res_disp_and_or_flag = 0

# if the res_disp_flag_and_or = 0 then it will be or if set to 1 it will be one
if res_disp_and_or_flag == 1:
    res_disp_and_or = 'or'
else:
    res_disp_and_or = 'and'
    
# setting tolerance convergence 
con_tol_res = 0.03 #0.01
con_tol_dis = 0.03

# to put on py_conect set flag to 1 here
py_proc_flag = 0

# only output the displacement for the model then set to 1
displacement_only_flag = 1

# turn of on creating procedure files:
procedure_flag = 0

# Beam element radius
beam_radius = fin_rad/5 # default fin_rad

# see beam
view_beam_flag = 1 #0
#*wire_element_display_style_normal

# flag to use RBE elements or nodal ties
RBEorTies_flag = 1 # 0 for RBE

# Selecting either Smooth-Sil 950 (option 0) or Moldstar 15 (option 1)
Sil_option = 1

bulk_auto_flag = 1 #(1 = use defined bulk values; 0 ignore and use auto)

# Option to automatically define died and retained node: (1 = on; 0 = off)
aut_node_flag = 1

Sil_Mooney_C10 = 260567.62
Sil_Mooney_C01 = 97549.81
Sil_Mooney_C20 = 57500.69
density = 1211
bulk = 679e3
if Sil_option ==1:
    Sil_Mooney_C10 = 36428.1
    Sil_Mooney_C01 = 18924.3
    Sil_Mooney_C20 = 3926.09
    density = 1139
    bulk = 119e3
	
Marc_script = """"""

Marc_script = Marc_script + """| Created by Marc Mentat 2021.4 (64bit)
*prog_option compatibility:prog_version:ment2021.4
*prog_analysis_class structural
*prog_use_current_job on
*set_default_length_unit meter
*set_model_length_unit meter
|
*theme dark
*set_save_formatted off *save_as_model \""""+root_dir+"""\\FinalOpt\\Marc\\Model1.mud\" yes
*change_directory "."
"""+root_dir+"""\\FinalOpt\\Marc
"""
#====================
py_proc = """"""

if py_proc_flag == 1:
    py_proc = """
*py_separate_process on
"""
#===================
    
Marc_script = Marc_script +py_proc+"""

*set_undo off
*set_update off
*model_navigator_update off
"""

# This is to select nodes to apply 
#==============================
#******************************

finger_accurate = """
"""
#******************************
#==============================

Marc_script = Marc_script + finger_accurate


Marc_script = Marc_script+ """
*prog_option import_nastran:show_log:off
"""

# Marc_script = Marc_script + """
# *import nastran \""""+root_dir+"""\\APEX\\Actuator.bdf\"\n
# *import nastran \""""+root_dir+"""\\APEX\\Surface.bdf\"\n
# """
    
Marc_script = Marc_script + """

*new_pre_defined_table linear_ramp_time
*table_name linear_ramp
*new_pre_defined_table linear_ramp_up_down_time
*table_name linear_ramp_up_down
*current_graphics_window model:1
*fill_view
*grid_u_min -0.02
*grid_u_max 0.02
*grid_u_spacing 0.02
*grid_u_min -0.1
*grid_u_max 0.1
*grid_v_min -0.1
*grid_v_max 0.1
*grid_v_spacing 0.02
*system_grid_display_on
@popup(set_control_pm,0)
@popdown(set_control_pm,0)
@main(geometric_properties)
*new_geometry *geometry_type mech_three_solid
*geometry_name geometry_actuator
*edit_geometry geometry_actuator
*add_geometry_elements
all_visible
# | End of List
*geometry_option assumedstrn:on
*new_mater standard *mater_option general:state:solid *mater_option general:skip_structural:off
*mater_option structural:type:mooney
*mater_name Mooney
*mater_param structural:mooney_c10 """+str(Sil_Mooney_C10)+"""
*mater_param structural:mooney_c01 """+str(Sil_Mooney_C01)+"""
*mater_param structural:mooney_c20 """+str(Sil_Mooney_C20)+"""
*add_mater_elements
Actuator
# | End of List
"""

Marc_script = Marc_script + """
*new_mater standard *mater_option general:state:solid *mater_option general:skip_structural:off
*mater_name Paper
*mater_param structural:youngs_modulus 0.3e9
*mater_param structural:poissons_ratio 0.33
*edit_mater Paper
*add_mater_elements
Paper
# | End of List

*new_cbody mesh *contact_option state:solid *contact_option skip_structural:off
*contact_body_name CB_Actuator
*add_contact_body_elements
Actuator
# | End of List

*new_cbody mesh *contact_option state:solid *contact_option skip_structural:off
*contact_body_name CB_Strain_limit
*add_contact_body_elements
Paper
# | End of List

*new_interact mesh:mesh *interact_option state_1:solid *interact_option state_2:solid
*interact_name Glued
*interact_option contact_type:glue
*interact_option delay_slide_off:on
*interact_option project_stress_free:on
*new_interact mesh:mesh *interact_option state_1:solid *interact_option state_2:solid
*interact_name Touch
*new_contact_table
*visible_set Actuator *redraw
*visible_set Strain_limit *redraw
*ctable_entry CB_Actuator CB_Actuator
*contact_table_option CB_Actuator CB_Actuator contact:on
*prog_string ctable:old_interact Touch *ctable_entry_interact CB_Actuator CB_Actuator

*new_apply *apply_type fixed_displacement
*apply_name Fixed_side
*apply_dof x *apply_dof_value x
*apply_dof y *apply_dof_value y
*apply_dof z *apply_dof_value z
*apply_dof rx *apply_dof_value rx
*apply_dof ry *apply_dof_value ry
*apply_dof rz *apply_dof_value rz
*select_method_user_box
*select_filter_none
*select_nodes """ +str(round(Origin_start-2e-3-0.5e-3,6))+ " " + str(round(Origin_start-2e-3+0.5e-3,6))+ " " + str(round(-(width/2)-0.5e-3,6))+ " " + str(round((width/2)+0.5e-3,6))+ " " + str(round(-5.500000000000e-03,6)) + " " + str(round(height + 0.5e-3,6)) + """
*add_apply_nodes
all_selected
# | End of List
*select_clear_nodes
*model_orientation_bottom
*fill_view
*auto_arrow_length 0.03
*redraw

"""
# *draw_clipped on
# *set_clip_plane_normal_comp y -1
# *set_clip_plane_normal_comp x 0

ch_sel_1 = ch_var1#3.31666666666666665#3.8
no_ch_max_1 = (L1-L0-ac_th_for-gap)/(2*ac_th_for+0.001+gap) # NB CORRECT
if ch_sel_1>no_ch_max_1:
    ch_sel_1 = no_ch_max_1

ch_width_1 = (L1-L0-ac_th_for-gap)/ch_sel_1 - 2*ac_th_for-gap 

rem_n = ch_sel_1-math.floor(ch_sel_1)

if rigid1 == 0:
    copy_text = "*select_method_user_box\n*select_filter_surface\n*select_faces "
    text_comb_1 = ''
    st_pt = L0
    for a in range(math.floor(ch_sel_1)):
        x_min = str(round(st_pt+(2*ac_th_for+ch_width_1+gap)*(a)+ac_th_for-meshing_size/2,7))
        x_max = str(round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(a)-ac_th_for+meshing_size/2,7))    
        y_min = str(-ch_width/2-meshing_size/2) 
        y_max = str(ch_width/2+meshing_size/2)
        z_min = str(-meshing_size/2)
        z_max = str(height-ac_th_up+meshing_size/2)
        text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
        temp = copy_text+text_temp
        text_comb_1 = text_comb_1 + temp 
    #print(text_comb_1)

    if rem_n ==0:
        vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + ac_th_for   ,7)
        varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + (2*ac_th_for+ch_width_1)   ,7)
        varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - ac_th_side-gap   ,7)
        vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7)
    else:
        vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + ac_th_for   ,7)
        varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + (2*ac_th_for+ch_width_1)    ,7)
        varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - ac_th_side - gap   ,7)
        vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1))   ,7)


        if varc <= vara:
            print('nothing')
        elif varc > vara:
            if varc > varb:
                x_min = str(round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for-meshing_size/2,7))
                                 #st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                x_max = str(round(varc+meshing_size/2,7))
                                 #st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for
                y_min = str(-ch_width/2-meshing_size/2) 
                y_max = str(ch_width/2+meshing_size/2)
                z_min = str(-meshing_size/2)
                z_max = str(height-ac_th_up+meshing_size/2)
                text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
                temp = copy_text+text_temp
                text_comb_1 = text_comb_1 + temp
                
            elif varc <= varb:
                x_min = str(round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for-meshing_size/2,7))
                                 #st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                x_max = str(round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for+ch_width_1+meshing_size/2,7))     
                y_min = str(-ch_width/2-meshing_size/2) 
                y_max = str(ch_width/2+meshing_size/2)
                z_min = str(-meshing_size/2)
                z_max = str(height-ac_th_up+meshing_size/2)
                text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
                temp = copy_text+text_temp
                text_comb_1 = text_comb_1 + temp 

        if varc <= vara:
            print('nothing')
        elif varc > vara:
            if varc <= vard:
                ca = st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                   #st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                cb = varc
                   # varc
                x_min = str(round(ca-meshing_size/2,6))
                x_max = str(round(cb+meshing_size/2,6))
                y_min = str(-ch_width/2-meshing_size/2) 
                y_max = str(ch_width/2+meshing_size/2)
                z_min = str(-meshing_size/2)
                z_max = str(height-ac_th_up+meshing_size/2)
                text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
                temp = copy_text+text_temp
                text_comb_1 = text_comb_1 + temp
            elif varc > vard:
                if varc <= varb:
                    ca = st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                        #st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                    cb = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for
                        #st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for
                    x_min = str(round(ca-meshing_size/2,6))
                    x_max = str(round(cb+meshing_size/2,6))
                    y_min = str(-ch_width/2-meshing_size/2) 
                    y_max = str(ch_width/2+meshing_size/2)
                    z_min = str(-meshing_size/2)
                    z_max = str(height-ac_th_up+meshing_size/2)
                    text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
                    temp = copy_text+text_temp
                    text_comb_1 = text_comb_1 + temp

                elif varc > varb:
                    ca = st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                        #st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                    cb = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for
                        #st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for
                    x_min = str(round(ca-meshing_size/2,6))
                    x_max = str(round(cb+meshing_size/2,6))
                    y_min = str(-ch_width/2-meshing_size/2) 
                    y_max = str(ch_width/2+meshing_size/2)
                    z_min = str(-meshing_size/2)
                    z_max = str(height-ac_th_up+meshing_size/2)
                    text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
                    temp = copy_text+text_temp
                    text_comb_1 = text_comb_1 + temp

    Marc_script = Marc_script + text_comb_1

Marc_script = Marc_script + """
*store_faces
MC_top
all_selected
*select_clear
"""


final_combined_script = """"""

text_final_mat = """"""
text_final_mat = text_final_mat + """
*edit_apply Ch1_Bot
*apply_dof_value p """+str(round(LHC_P*300000,2))+"""
*edit_rbe2 rbe2_1 *remove_current_rbe2
*edit_rbe2 rbe2_2 *remove_current_rbe2
*edit_rbe2 rbe2_3 *remove_current_rbe2
*remove_elements
all_existing
*remove_nodes
all_existing

*change_directory .
D:\\SB_SIMULATIONS\\FinalResults\\Chamb3_TwoLat\\"""+str('S')+"""_"""+ str(dim_c) +"""
*set_save_formatted off *save_as_model "..\\"""+ str('S')+"""_"""+str(dim_c)+"""\\"""+ str(dim)+"""_"""+str(dim_c) +""".mud" yes
*prog_option import_nastran:show_log:off
*select_mode_and
*select_method_single
*select_filter_none
*import nastran \"C:\\Users\\18507522\\FinalSBSpring\\CODE\\3DModelSave\\3BlockActuatorFinal\\"""+ str('Sel')+str(dim_c) +"""_2D.bdf\"
*select_elements_geometry Elements_and_Element_Properties_for:_Paper
*store_elements Paper
all_selected
*select_clear_elements
*edit_geometry Elements_and_Element_Properties_for:_Paper
*edit_geometry Elements_and_Element_Properties_for:_Paper *remove_current_geometry
*import nastran \"C:\\Users\\18507522\\FinalSBSpring\\CODE\\3DModelSave\\3BlockActuatorFinal\\"""+ str('Sel')+str(dim_c) +"""_3D.bdf\"
*select_elements_geometry Elements_and_Element_Properties_for:_Actuator
*store_elements Actuator
all_selected
*select_clear_elements
*edit_geometry Elements_and_Element_Properties_for:_Actuator
*edit_geometry Elements_and_Element_Properties_for:_Actuator *remove_current_geometry

"""+text_comb_1+"""
*edit_apply Ch1_Top
*add_apply_faces
all_selected
*select_clear_faces
"""

#print('*edit_apply Ch1_Bot')
#print('*apply_dof_value p ',round(LHC_P*300000,2))
#print('*edit_rbe2 rbe2_1 *remove_current_rbe2')
#print('*edit_rbe2 rbe2_2 *remove_current_rbe2')
#print('*edit_rbe2 rbe2_3 *remove_current_rbe2')
#print('*remove_elements')
#print('all_existing')
#print('*remove_nodes')
#print('all_existing')
#print('theme dark')
code_up = """
*prog_option import_nastran:show_log:off
*select_mode_and
*select_method_single
*select_filter_none
*import nastran \"C:\\Users\\18507522\\FinalSBSpring\\CODE\\3DModelSave\\3BlockActuator\\Sel"""+ str(samp) +"""_2D.bdf\"
*select_elements_geometry Elements_and_Element_Properties_for:_Paper
*store_elements Paper
all_selected
*select_clear_elements
*edit_geometry Elements_and_Element_Properties_for:_Paper
*edit_geometry Elements_and_Element_Properties_for:_Paper *remove_current_geometry
*import nastran \"C:\\Users\\18507522\\FinalSBSpring\\CODE\\3DModelSave\\3BlockActuator\\Sel"""+ str(samp) +"""_3D.bdf\"
*select_elements_geometry Elements_and_Element_Properties_for:_Actuator
*store_elements Actuator
all_selected
*select_clear_elements
*edit_geometry Elements_and_Element_Properties_for:_Actuator
*edit_geometry Elements_and_Element_Properties_for:_Actuator *remove_current_geometry
"""
#print(code_up)
#print('')
#print(text_comb_1)
#print("*edit_apply Ch1_Top")
#print("*add_apply_faces")
#print("all_selected")
#print("*select_clear_faces")

ch_sel_1 = ch_var1#3.31666666666666665#3.8
no_ch_max_1 = (L1-L0-ac_th_for-gap)/(2*ac_th_for+0.001+gap) # NB CORRECT
if ch_sel_1>no_ch_max_1:
    ch_sel_1 = no_ch_max_1

ch_width_1 = (L1-L0-ac_th_for-gap)/ch_sel_1 - 2*ac_th_for-gap 

rem_n = ch_sel_1-math.floor(ch_sel_1)

if rigid1 == 0:
    copy_text = "*select_method_user_box\n*select_filter_surface\n*select_faces "
    text_comb_1 = ''
    st_pt = L0
    for a in range(math.floor(ch_sel_1)):
        x_min = str(round(st_pt+(2*ac_th_for+ch_width_1+gap)*(a)+ac_th_for-meshing_size/2,7))
        x_max = str(round(st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(a)-ac_th_for+meshing_size/2,7))    
        y_min = str(-ch_width/2-meshing_size/2) 
        y_max = str(ch_width/2+meshing_size/2)
        z_min = str(-0.004-(height-ac_th_up+meshing_size/2))
        z_max = str(meshing_size/2-0.004)
        text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
        temp = copy_text+text_temp
        text_comb_1 = text_comb_1 + temp 
    #print(text_comb_1)

    if rem_n ==0:
        vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + ac_th_for   ,7)
        varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1-1)) + (2*ac_th_for+ch_width_1)   ,7)
        varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) - ac_th_side-gap   ,7)
        vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))   ,7)
    else:
        vara = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + ac_th_for   ,7)
        varb = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1)) + (2*ac_th_for+ch_width_1)    ,7)
        varc = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - (2*ac_th_for+ch_width_1+gap)*(1-rem_n)   ,7)
        vard = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1)) - ac_th_side - gap   ,7)
        vare = round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1+1))   ,7)


        if varc <= vara:
            print('nothing')
        elif varc > vara:
            if varc > varb:
                x_min = str(round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for-meshing_size/2,7))
                                 #st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                x_max = str(round(varc+meshing_size/2,7))
                                 #st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for
                y_min = str(-ch_width/2-meshing_size/2) 
                y_max = str(ch_width/2+meshing_size/2)
                z_min = str(-0.004-(height-ac_th_up+meshing_size/2))
                z_max = str(meshing_size/2-0.004)
                text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
                temp = copy_text+text_temp
                text_comb_1 = text_comb_1 + temp
                
            elif varc <= varb:
                x_min = str(round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for-meshing_size/2,7))
                                 #st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                x_max = str(round(st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for+ch_width_1+meshing_size/2,7))     
                y_min = str(-ch_width/2-meshing_size/2) 
                y_max = str(ch_width/2+meshing_size/2)
                z_min = str(-0.004-(height-ac_th_up+meshing_size/2))
                z_max = str(meshing_size/2-0.004)
                text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
                temp = copy_text+text_temp
                text_comb_1 = text_comb_1 + temp 

        if varc <= vara:
            print('nothing')
        elif varc > vara:
            if varc <= vard:
                ca = st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                   #st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                cb = varc
                   # varc
                x_min = str(round(ca-meshing_size/2,6))
                x_max = str(round(cb+meshing_size/2,6))
                y_min = str(-ch_width/2-meshing_size/2) 
                y_max = str(ch_width/2+meshing_size/2)
                z_min = str(-0.004-(height-ac_th_up+meshing_size/2))
                z_max = str(meshing_size/2-0.004)
                text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
                temp = copy_text+text_temp
                text_comb_1 = text_comb_1 + temp
            elif varc > vard:
                if varc <= varb:
                    ca = st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                        #st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                    cb = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for
                        #st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for
                    x_min = str(round(ca-meshing_size/2,6))
                    x_max = str(round(cb+meshing_size/2,6))
                    y_min = str(-ch_width/2-meshing_size/2) 
                    y_max = str(ch_width/2+meshing_size/2)
                    z_min = str(-0.004-(height-ac_th_up+meshing_size/2))
                    z_max = str(meshing_size/2-0.004)
                    text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
                    temp = copy_text+text_temp
                    text_comb_1 = text_comb_1 + temp

                elif varc > varb:
                    ca = st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                        #st_pt+(2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))+ac_th_for
                    cb = st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for
                        #st_pt+(2*ac_th_for+ch_width_1) + (2*ac_th_for+ch_width_1+gap)*(math.floor(ch_sel_1))-ac_th_for
                    x_min = str(round(ca-meshing_size/2,6))
                    x_max = str(round(cb+meshing_size/2,6))
                    y_min = str(-ch_width/2-meshing_size/2) 
                    y_max = str(ch_width/2+meshing_size/2)
                    z_min = str(-0.004-(height-ac_th_up+meshing_size/2))
                    z_max = str(meshing_size/2-0.004)
                    text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
                    temp = copy_text+text_temp
                    text_comb_1 = text_comb_1 + temp

    Marc_script = Marc_script + text_comb_1

Marc_script = Marc_script + """
*store_faces
MC_bot
all_selected
*select_clear
"""

text_final_mat = text_final_mat + """
"""+text_comb_1+"""
*edit_apply Ch1_Bot
*add_apply_faces
all_selected
*select_clear_faces

*edit_geometry geom1
*add_geometry_elements
Paper
*edit_mater Mat_Paper
*add_mater_elements
Paper
*edit_mater Mooney
*add_mater_elements
Actuator
*edit_contact_body Con_Actuator
*add_contact_body_elements
Actuator

*element_type 84
Actuator
*element_type 139
Paper

*edit_job job 
*job_option assem_recov_multi_threading:on
*job_param assem_recov_nthreads 4
*job_option nsolver_procs_serial:on
*job_param nthreads 4

*new_rbe2
*rbe2_tied_dof x
*rbe2_tied_dof y
*rbe2_tied_dof z
*new_rbe2
*rbe2_tied_dof x
*rbe2_tied_dof y
*rbe2_tied_dof z
*new_rbe2
*rbe2_tied_dof x
*rbe2_tied_dof y
*rbe2_tied_dof z
*rbe2_tied_dof rx
*rbe2_tied_dof rz
*select_clear
*select_reset
"""

final_combined_script = final_combined_script + text_final_mat
#print(text_comb_1)
#print("*edit_apply Ch1_Bot")
#print("*add_apply_faces")
#print("all_selected")
#print("*select_clear_faces")
#print(text_final_mat)

# correct 1
final_adapt =""""""
temp = """
*select_reset
*select_method_user_box
*select_filter_surface
"""
sel_side = """"""
st_pt = L0
# width height bredth
x_min = str(round(st_pt+2*ac_th_for+ch_width_1-2*meshing_size-meshing_size/2,7))
x_max = str(round(st_pt+2*ac_th_for+ch_width_1+gap+2*meshing_size+meshing_size/2,7))    
y_min = str(-width/2-meshing_size/2) 
y_max = str(width/2+meshing_size/2)
z_min = str(-0.004-(height+meshing_size/2))
z_max = str(3*meshing_size-0.004-(height+meshing_size/2))

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"*select_elements "+text_temp

#print(temp)
final_adapt = final_adapt+temp

# correct 2
temp = """
*select_reset
*select_method_user_box
*select_filter_surface
"""
sel_side = """"""
st_pt = L0
# width height bredth
x_min = str(round(st_pt+2*ac_th_for+ch_width_1-2*meshing_size-meshing_size/2,7))
x_max = str(round(st_pt+2*ac_th_for+ch_width_1+gap+2*meshing_size+meshing_size/2,7))    
y_min = str(-width/2-meshing_size/2) 
y_max = str(width/2+meshing_size/2)
z_min = str(height-meshing_size/2-meshing_size*3)
z_max = str(height+meshing_size/2)

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"*select_elements "+text_temp
#print(temp)

final_adapt = final_adapt+temp

# correct 3
temp = """
*select_reset
*select_method_user_box
*select_filter_surface
"""
sel_side = """"""
st_pt = L0
# width height bredth
x_min = str(round(gap+2*ac_th_for+ch_width_1+st_pt+2*ac_th_for+ch_width_1-2*meshing_size-meshing_size/2,7))
x_max = str(round(+gap+2*ac_th_for+ch_width_1+st_pt+2*ac_th_for+ch_width_1+gap+2*meshing_size+meshing_size/2,7))    
y_min = str(-width/2-meshing_size/2) 
y_max = str(width/2+meshing_size/2)
z_min = str(-0.004-(height+meshing_size/2))
z_max = str(3*meshing_size-0.004-(height+meshing_size/2))

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"*select_elements "+text_temp

#print(temp)

final_adapt = final_adapt+temp

# correct 4
temp = """
*select_reset
*select_method_user_box
*select_filter_surface
"""
sel_side = """"""
st_pt = L0
# width height bredth
x_min = str(round(gap+2*ac_th_for+ch_width_1+st_pt+2*ac_th_for+ch_width_1-2*meshing_size-meshing_size/2,7))
x_max = str(round(+gap+2*ac_th_for+ch_width_1+st_pt+2*ac_th_for+ch_width_1+gap+2*meshing_size+meshing_size/2,7))     
y_min = str(-width/2-meshing_size/2) 
y_max = str(width/2+meshing_size/2)
z_min = str(height-meshing_size/2-meshing_size*3)
z_max = str(height+meshing_size/2)

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"*select_elements "+text_temp
#print(temp)

final_adapt = final_adapt+temp

temp = """
*select_reset
*select_method_user_box
*select_filter_surface
"""
sel_side = """"""
st_pt = L0
# width height bredth
x_min = str(round(gap/2+st_pt+2*ac_th_for+ch_width_1-2*meshing_size,7))
x_max = str(round(gap/2+st_pt+2*ac_th_for+ch_width_1+2*meshing_size,7))    
y_min = str(-width/2-meshing_size/2) 
y_max = str(width/2+meshing_size/2)
z_min = str(-meshing_size*2-0.004-0.003-meshing_size/2)
z_max = str(meshing_size*2-0.004-0.003+meshing_size/2)# 0.003 is gap height

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"*select_elements "+text_temp

#print(temp)

final_adapt = final_adapt+temp

temp = """
*select_reset
*select_method_user_box
*select_filter_surface
"""
sel_side = """"""
st_pt = L0
# width height bredth
x_min = str(round(gap/2+st_pt+2*ac_th_for+ch_width_1-2*meshing_size,7))
x_max = str(round(gap/2+st_pt+2*ac_th_for+ch_width_1+2*meshing_size,7))    
y_min = str(-width/2-meshing_size/2) 
y_max = str(width/2+meshing_size/2)
z_min = str(-meshing_size*2+0.003-meshing_size/2)
z_max = str(meshing_size*2+0.003+meshing_size/2)

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"*select_elements "+text_temp

#print(temp)

final_adapt = final_adapt+temp

temp = """
*select_reset
*select_method_user_box
*select_filter_surface
"""
sel_side = """"""
st_pt = L0
# width height bredth
x_min = str(round(gap*1.5+st_pt+4*ac_th_for+2*ch_width_1-2*meshing_size,7))
x_max = str(round(gap*1.5+st_pt+4*ac_th_for+2*ch_width_1+2*meshing_size,7))    
y_min = str(-width/2-meshing_size/2) 
y_max = str(width/2+meshing_size/2)
z_min = str(-meshing_size*2-0.004-0.003-meshing_size/2)
z_max = str(meshing_size*2-0.004-0.003+meshing_size/2)# 0.003 is gap height

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"*select_elements "+text_temp

#print(temp)

final_adapt = final_adapt+temp

temp = """
*select_reset
*select_method_user_box
*select_filter_surface
"""
sel_side = """"""
st_pt = L0
# width height bredth
x_min = str(round(gap*1.5+st_pt+4*ac_th_for+2*ch_width_1-2*meshing_size,7))
x_max = str(round(gap*1.5+st_pt+4*ac_th_for+2*ch_width_1+2*meshing_size,7))    
y_min = str(-width/2-meshing_size/2) 
y_max = str(width/2+meshing_size/2)
z_min = str(-meshing_size*2+0.003-meshing_size/2)
z_max = str(meshing_size*2+0.003+meshing_size/2)

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"*select_elements "+text_temp

#print(temp)

final_adapt = final_adapt+temp

final_adapt = final_adapt+"""
*edit_adapt adapt1
*add_adapt_elements
all_selected
*select_clear_elements
"""
final_combined_script = final_combined_script + final_adapt

temp = """
*select_reset
*select_clear
*select_method_user_box
*select_filter_surface
"""
sel_side = """"""
st_pt = L0
# width height bredth
x_min = str(round(st_pt-meshing_size/2,7))
x_max = str(round(st_pt+meshing_size/2,7))    
y_min = str(-width/2-meshing_size/2) 
y_max = str(width/2+meshing_size/2)
z_min = str(-0.004-(height+meshing_size/2))
z_max = str(height+meshing_size/2)

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"*select_nodes "+text_temp
#text_comb_1 = text_comb_1 + temp 
#print(text_comb_1)


temp = temp+ """
*edit_apply Fixed_side
*add_apply_nodes
all_selected

*select_reset
*select_clear
*select_method_user_box
*select_filter_surface
"""
text_temp = """"""
x_min = str(round(st_pt-meshing_size/2,7))
x_max = str(round(L1+meshing_size/2,7))    
y_min = str(0-meshing_size/2) 
y_max = str(0+meshing_size/2)
z_min = str(-0.004-(height+meshing_size/2))
z_max = str(height+meshing_size/2)

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"*select_nodes "+text_temp

temp = temp +"""
*edit_apply symmetry_side
*add_apply_nodes
all_selected

*select_reset
*select_clear
"""


temp = temp+ """
*edit_apply Fixed_side
*add_apply_nodes
all_selected

*select_reset
*select_clear
*select_method_user_box
*select_filter_surface
"""
text_temp = """"""
x_min = str(round(L1-1e-3-meshing_size/2,7))
x_max = str(round(L1-1e-3+meshing_size/2,7))    
y_min = str(-width/2-meshing_size/2) 
y_max = str(width/2+meshing_size/2)
z_min = str(-0.004-(height+meshing_size/2))
z_max = str(height+meshing_size/2)

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"*select_nodes "+text_temp


text_temp = """"""
x_min = str(round(L1-1e-3-meshing_size/2,7))
x_max = str(round(L1-1e-3+meshing_size/2,7))    
y_min = str(-meshing_size/2) 
y_max = str(meshing_size/2)
z_min = str(-0.002-(meshing_size/2))
z_max = str(-0.002+meshing_size/2)

text_temp = x_min + ' '+ x_max + ' ' + y_min + ' '+y_max + ' '+z_min + ' '+z_max + '\n'
temp = temp+"""
*select_mode_except
"""+"*select_nodes "+text_temp + """
*select_reset



@popup(meshgen_pm,0)
"""

final_combined_script = final_combined_script + temp
#print(temp)
temp = """"""

#print(final_combined_script)
# # NOTE NB THIS PRINTS THE FINAL CODE TILL 
# text_hold_fin
with open(root_dir+"\\"+"ThreeActuator.proc", 'w') as f:
    f.write(final_combined_script)
    f.close()

print(root_dir+"\\"+"ThreeActuator.proc")