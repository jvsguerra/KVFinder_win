####################################################################################################
#    This file is part of KVFinder.                                                                #
#                                                                                                  #
#    KVFinder is free software: you can redistribute it and/or modify                              #
#    it under the terms of the GNU General Public License as published by                          #
#    the Free Software Foundation, either version 3 of the License, or                             #
#    (at your option) any later version.                                                           #
#                                                                                                  #
#    KVFinder is distributed in the hope that it will be useful,                                   #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                                #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 #
#    GNU General Public License for more details.                                                  #
#                                                                                                  #
#    You should have received a copy of the GNU General Public License                             #
#    along with KVFinder.  If not, see <http://www.gnu.org/licenses/>.                             #
#                                                                                                  #
#                                                                                                  #
#    The KVFinder software was developed by:                                                       #
#    Saulo Henrique Pires de Oliveira                                                              #
#    Felipe Augusto Nunes Ferraz                                                                   #
#    Rodrigo Vargas Honorato                                                                       #
#    Jose Xavier Neto                                                                              #
#    Tiago Jose Paschoal Sobreira                                                                  #
#    Paulo Sergio Lopes de Oliveira                                                                #
#                                                                                                  #
#    National Center of Energy and Material Research - CNPEM                                       #
#    National Laboratory of Biosciences - LNBio                                                    #
#    Campinas, Brazil - P.O. Box 6192 - CEP 13083-970, Campinas - SP                               #
#                                                                                                  #
#    Contact: paulo.oliveira@lnbio.cnpem.br                                                        #
#    KVFinder website: http://lnbio.cnpem.br/bioinformatics/main/software                          #
####################################################################################################


####################################################################################################
#This is the source code of the Windows KVFinder Plugin Interface for Pymol. It was developed using#
#TK and Python.                                                                                    #
#                                                                                                  #
#                                                                                                  #
#After the installation of KVFinder, a copy of this file will be in the pymol directory.           #
#Changes in this file are not advised, as it controls all of the KVFinder features.                #
#                                                                                                  #
#                                                                                                  #
####################################################################################################

import tkSimpleDialog
import tkMessageBox
from Tkinter import *
import math
from pymol import cmd
from pymol.cgo import *
import pymol
import sys, zlib, commands
import Pmw
import tkFileDialog


#Initialization of all the global variables in the default values

box_name = "box"
x = 0.0
y = 0.0
z = 0.0
lim1 = 0.0
lim2 = 0.0
lim3 = 0.0
lim4 = 0.0
lim5 = 0.0
lim6 = 0.0
angle1 = 0.0
angle2 = 0.0

ligandAd = BooleanVar()
ligandAd.set(0)
StepRed = BooleanVar()
StepRed.set(1)
WholeProt = BooleanVar()
WholeProt.set(0)
OverW = BooleanVar()
OverW.set(0)
ChargeDisplay = BooleanVar()
ChargeDisplay.set(0)
BoxMode = BooleanVar()
BoxMode.set(1)


infile = StringVar()
path = StringVar()
#path.set(os.getcwd())

#if os.path.isdir(path.get()) == 0:
#    os.makedirs(path.get())

#os.chdir(path.get())
path.set(os.path.expanduser("~"))
if os.path.isdir(path.get()+"\\KVFinder") == 0:
    os.makedirs(path.get()+'\\KVFinder')
path.set(path.get()+'\\KVFinder\\')
os.chdir(path.get())
outfile = StringVar()
outfileBasic = StringVar()
outfile.set('out')
outfileBasic.set('out')
last_object = StringVar()
last_object.set("")
last_input = StringVar()
last_input.set("")
residues = {}


# this creates a menu entry
# and calls a function to execute the commands
def __init__(self):
   self.menuBar.addmenuitem('Plugin', 'command',
                            'KVFinder',
                            label = 'KVFinder',
                            command = lambda s=self : testDialog())


#Save the Parameters.txt file, used in the kvfinder
def saveBox(d, s, p, v, surf, lf, dc, o, p_out,par):

        f_out = open("Parameters.txt", "w")
        f_out.write(d.get()+"\n"+infile.get()+"\n"+o.get()+".KVFinder.output\n"+str(WholeProt.get())+"\n"+s.get()+"\n"+str(StepRed.get())+"\n\n")

        z1 = lim1*math.sin(angle2) + (lim3)*math.sin(angle1)*math.cos(angle2) - (lim5)*math.cos(angle1)*math.cos(angle2) + z
        y1 = (-lim3)*math.cos(angle1) + (-lim5)*math.sin(angle1) + y
        x1 = (-lim1)*math.cos(angle2) - (-lim3)*math.sin(angle1)*math.sin(angle2) + (-lim5)*math.cos(angle1)*math.sin(angle2) + x

        f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

        z1 = -(lim2)*math.sin(angle2) - (-lim3)*math.sin(angle1)*math.cos(angle2) + (-lim5)*math.cos(angle1)*math.cos(angle2) + z
        y1 = (-lim3)*math.cos(angle1) + (-lim5)*math.sin(angle1) + y
        x1 = (lim2)*math.cos(angle2) - (-lim3)*math.sin(angle1)*math.sin(angle2) + (-lim5)*math.cos(angle1)*math.sin(angle2) + x

        f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

        z1 = -(-lim1)*math.sin(angle2) - (lim4)*math.sin(angle1)*math.cos(angle2) + (-lim5)*math.cos(angle1)*math.cos(angle2) + z
        y1 = (lim4)*math.cos(angle1) + (-lim5)*math.sin(angle1) + y
        x1 = (-lim1)*math.cos(angle2) - (lim4)*math.sin(angle1)*math.sin(angle2) + (-lim5)*math.cos(angle1)*math.sin(angle2) + x

        f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

        z1 = -(-lim1)*math.sin(angle2) - (-lim3)*math.sin(angle1)*math.cos(angle2) + (lim6)*math.cos(angle1)*math.cos(angle2) + z
        y1 = (-lim3)*math.cos(angle1) + (lim6)*math.sin(angle1) + y
        x1 = (-lim1)*math.cos(angle2) - (-lim3)*math.sin(angle1)*math.sin(angle2) + (lim6)*math.cos(angle1)*math.sin(angle2) + x

        f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

        z1 = (lim1+float(p_out.get()))*math.sin(angle2) + (lim3+float(p_out.get()))*math.sin(angle1)*math.cos(angle2) - (lim5+float(p_out.get()))*math.cos(angle1)*math.cos(angle2) + z
        y1 = (-lim3-float(p_out.get()))*math.cos(angle1) + (-lim5-float(p_out.get()))*math.sin(angle1) + y
        x1 = (-lim1-float(p_out.get()))*math.cos(angle2) - (-lim3-float(p_out.get()))*math.sin(angle1)*math.sin(angle2) + (-lim5-float(p_out.get()))*math.cos(angle1)*math.sin(angle2) + x

        f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

        z1 = -(lim2+float(p_out.get()))*math.sin(angle2) - (-lim3-float(p_out.get()))*math.sin(angle1)*math.cos(angle2) + (-lim5-float(p_out.get()))*math.cos(angle1)*math.cos(angle2) + z
        y1 = (-lim3-float(p_out.get()))*math.cos(angle1) + (-lim5-float(p_out.get()))*math.sin(angle1) + y
        x1 = (lim2+float(p_out.get()))*math.cos(angle2) - (-lim3-float(p_out.get()))*math.sin(angle1)*math.sin(angle2) + (-lim5-float(p_out.get()))*math.cos(angle1)*math.sin(angle2) + x

        f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

        z1 = -(-lim1-float(p_out.get()))*math.sin(angle2) - (lim4+float(p_out.get()))*math.sin(angle1)*math.cos(angle2) + (-lim5-float(p_out.get()))*math.cos(angle1)*math.cos(angle2) + z
        y1 = (lim4+float(p_out.get()))*math.cos(angle1) + (-lim5-float(p_out.get()))*math.sin(angle1) + y
        x1 = (-lim1-float(p_out.get()))*math.cos(angle2) - (lim4+float(p_out.get()))*math.sin(angle1)*math.sin(angle2) + (-lim5-float(p_out.get()))*math.cos(angle1)*math.sin(angle2) + x

        f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

        z1 = -(-lim1-float(p_out.get()))*math.sin(angle2) - (-lim3-float(p_out.get()))*math.sin(angle1)*math.cos(angle2) + (lim6+float(p_out.get()))*math.cos(angle1)*math.cos(angle2) + z
        y1 = (-lim3-float(p_out.get()))*math.cos(angle1) + (lim6+float(p_out.get()))*math.sin(angle1) + y
        x1 = (-lim1-float(p_out.get()))*math.cos(angle2) - (-lim3-float(p_out.get()))*math.sin(angle1)*math.sin(angle2) + (lim6+float(p_out.get()))*math.cos(angle1)*math.sin(angle2) + x

        f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

        f_out.write("\n"+p.get()+"\n"+p_out.get()+"\n"+str(BoxMode.get())+"\n"+v.get()+"\n"+str(surf.get())+"\n0\n"+str(ligandAd.get())+"\n"+lf.get()+"\n"+dc.get())

        f_out.close()

        lf_2 = Listbox(par, exportselection=0)
        lf_2.insert(0,lf.get())
        lf_2.selection_set(0)

        saveNewBox(d, s, p, v, surf, lf_2, dc, o, p_out,par)

#Save the NewParameters.txt file, used in the kvfinder
def saveNewBox(d, s, p, v, surf, lf, dc, o, p_out, par):

    global x, y, z

    if( lim1==0 and lim2==0 and lim3==0 and lim4==0 and lim5==0 and lim6==0 and WholeProt.get()==0):
        tkMessageBox.showerror("Warning", "Search space not defined. Please, define a search space by selecting some residues and clicking in Draw New Box or use the Whole Protein mode", parent=par)
        return

    if not s.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid step size", parent=par)
        return

    if not p.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid probe size", parent=par)
        return

    if not p_out.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid probe out size", parent=par)
        return

    if (WholeProt.get() == 1 or BoxMode.get() == 1) and StepRed.get() == 0 and float(s.get()) <= 1 and  float(p_out.get()) >= 20:
        if not tkMessageBox.askyesno("Warning", "Given the Probe Out defined, this operation may take a really long time. Should it proceed?", parent=par):
            return

    if not v.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid minimum cavity volume", parent=par)
        return

    if not dc.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid distance cutoff", parent=par)
        return


    if not os.path.isfile(d.get()):
        tkMessageBox.showerror("Warning", "Invalid Dictionary Filename", parent=par)
        return


    if ligandAd.get() == 1 and not lf.curselection():
        tkMessageBox.showerror("Warning", "Ligand Adjustment selected with invalid ligand filename", parent=par)
        return

    control = 1
    if ligandAd.get() == 1:
        for ob in cmd.get_names("all"):
            if ob == lf.get(lf.curselection()):
                control = 0

        if os.path.isfile(lf.get(lf.curselection())) == 1:
            control = 0

        if control == 1:
            tkMessageBox.showerror("Warning", "Ligand Adjustment selected with invalid ligand", parent=par)
            return

    if ligandAd.get() == 1 and float(dc.get()) <= 0:
        tkMessageBox.showerror("Warning", "Ligand Adjustment selected with invalid Distance Cutoff", parent=par)
        return

    if WholeProt.get() == 1 and (float(p.get()) >= float(p_out.get())):
        tkMessageBox.showerror("Warning", "Invalid probe sizes for Whole protein mode. Probe Out must be bigger than probe", parent=par)
        return

    if BoxMode.get() == 1 and (float(p.get()) >= float(p_out.get())):
        tkMessageBox.showerror("Warning", "Invalid probe sizes for Probe Out Adjustment mode. Probe Out must be bigger than probe", parent=par)
        return


    if infile.get() == "":
        infile.set("-")
    if outfile.get() == "":
        outfile.set("-")

    f_out = open("KVParameters.txt", "w")
    f_out.write("#This is a configuration file used by KVFinder. The configurable parameters are identified by >. All parameters must have a value.\n\n")
    f_out.write("#Path for the KVFinder dictionary, with atoms information. For a example, check the default dictionary file on the KVFinder folder.\n")
    f_out.write(">DIC_PATH\n"+d.get()+"\n\n#Path for the input PDB. When using KVFinder on the command line mode, it must be specified using a command line parameter.\n"
    +">INPUT_PATH\n"+infile.get()+"\n\n#Path and base name for the output files (A PDB and results file). On the command line mode, it will be defined automatically.\n"
    +">OUTPUT_PATH\n"+o.get()+".KVFinder.output\n\n#Whole Protein mode, work as boolean (1 on and 0 off). Defines the search space as the whole protein.\n"
    +">WHOLE_MODE\n"+str(WholeProt.get())+"\n\n#Defines the size between grid points. Directly affects the precision. Also has a effect on the running time.\n"
    +">GRID\n"+s.get()+"\n\n#If set (1), it limits the number of grid points, avoiding extremely long analysis. But may be set off (0) for detailed analysis.\n"
    +">STEP_RED\n"+str(StepRed.get())+"\n\n>BOX_COORDINATES\n#Coordinates of the box vertices, that defines the search space when the whole protein mode is off.\n#The first 4 points represents the visible box on the interface, while the rest are used for internal computations only.\n#When using the command line mode, this values can be easily defined using the pymol interface and the generate parameters option.\n")

    z1 = lim1*math.sin(angle2) + (lim3)*math.sin(angle1)*math.cos(angle2) - (lim5)*math.cos(angle1)*math.cos(angle2) + z
    y1 = (-lim3)*math.cos(angle1) + (-lim5)*math.sin(angle1) + y
    x1 = (-lim1)*math.cos(angle2) - (-lim3)*math.sin(angle1)*math.sin(angle2) + (-lim5)*math.cos(angle1)*math.sin(angle2) + x

    f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

    z1 = -(lim2)*math.sin(angle2) - (-lim3)*math.sin(angle1)*math.cos(angle2) + (-lim5)*math.cos(angle1)*math.cos(angle2) + z
    y1 = (-lim3)*math.cos(angle1) + (-lim5)*math.sin(angle1) + y
    x1 = (lim2)*math.cos(angle2) - (-lim3)*math.sin(angle1)*math.sin(angle2) + (-lim5)*math.cos(angle1)*math.sin(angle2) + x

    f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

    z1 = -(-lim1)*math.sin(angle2) - (lim4)*math.sin(angle1)*math.cos(angle2) + (-lim5)*math.cos(angle1)*math.cos(angle2) + z
    y1 = (lim4)*math.cos(angle1) + (-lim5)*math.sin(angle1) + y
    x1 = (-lim1)*math.cos(angle2) - (lim4)*math.sin(angle1)*math.sin(angle2) + (-lim5)*math.cos(angle1)*math.sin(angle2) + x

    f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

    z1 = -(-lim1)*math.sin(angle2) - (-lim3)*math.sin(angle1)*math.cos(angle2) + (lim6)*math.cos(angle1)*math.cos(angle2) + z
    y1 = (-lim3)*math.cos(angle1) + (lim6)*math.sin(angle1) + y
    x1 = (-lim1)*math.cos(angle2) - (-lim3)*math.sin(angle1)*math.sin(angle2) + (lim6)*math.cos(angle1)*math.sin(angle2) + x

    f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

    z1 = (lim1+float(p_out.get()))*math.sin(angle2) + (lim3+float(p_out.get()))*math.sin(angle1)*math.cos(angle2) - (lim5+float(p_out.get()))*math.cos(angle1)*math.cos(angle2) + z
    y1 = (-lim3-float(p_out.get()))*math.cos(angle1) + (-lim5-float(p_out.get()))*math.sin(angle1) + y
    x1 = (-lim1-float(p_out.get()))*math.cos(angle2) - (-lim3-float(p_out.get()))*math.sin(angle1)*math.sin(angle2) + (-lim5-float(p_out.get()))*math.cos(angle1)*math.sin(angle2) + x

    f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

    z1 = -(lim2+float(p_out.get()))*math.sin(angle2) - (-lim3-float(p_out.get()))*math.sin(angle1)*math.cos(angle2) + (-lim5-float(p_out.get()))*math.cos(angle1)*math.cos(angle2) + z
    y1 = (-lim3-float(p_out.get()))*math.cos(angle1) + (-lim5-float(p_out.get()))*math.sin(angle1) + y
    x1 = (lim2+float(p_out.get()))*math.cos(angle2) - (-lim3-float(p_out.get()))*math.sin(angle1)*math.sin(angle2) + (-lim5-float(p_out.get()))*math.cos(angle1)*math.sin(angle2) + x

    f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

    z1 = -(-lim1-float(p_out.get()))*math.sin(angle2) - (lim4+float(p_out.get()))*math.sin(angle1)*math.cos(angle2) + (-lim5-float(p_out.get()))*math.cos(angle1)*math.cos(angle2) + z
    y1 = (lim4+float(p_out.get()))*math.cos(angle1) + (-lim5-float(p_out.get()))*math.sin(angle1) + y
    x1 = (-lim1-float(p_out.get()))*math.cos(angle2) - (lim4+float(p_out.get()))*math.sin(angle1)*math.sin(angle2) + (-lim5-float(p_out.get()))*math.cos(angle1)*math.sin(angle2) + x

    f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

    z1 = -(-lim1-float(p_out.get()))*math.sin(angle2) - (-lim3-float(p_out.get()))*math.sin(angle1)*math.cos(angle2) + (lim6+float(p_out.get()))*math.cos(angle1)*math.cos(angle2) + z
    y1 = (-lim3-float(p_out.get()))*math.cos(angle1) + (lim6+float(p_out.get()))*math.sin(angle1) + y
    x1 = (-lim1-float(p_out.get()))*math.cos(angle2) - (-lim3-float(p_out.get()))*math.sin(angle1)*math.sin(angle2) + (lim6+float(p_out.get()))*math.cos(angle1)*math.sin(angle2) + x

    f_out.write(str(x1)+" "+str(y1)+" "+str(z1)+"\n\n")

    f_out.write("#KVFinder works with a two sized probe system. A smaller probe, the Probe In, and a bigger one, the Probe Out, rolls around the protein. The points reached by the Probe In but not by the Probe Out are considered cavity points.\n#The Probe In radius.\n>PROBE_IN\n"+p.get()+"\n\n#The Probe Out radius\n>PROBE_OUT\n"
    +p_out.get()+"\n\n#This boolean defines if the Probe Out will be used (1) or not (0). Using only one probe might be useful for fast evaluations of areas. Volume results will be directly affected by the search space defined.\n>PROBE_OUT_ADJUSTMENT\n"+str(BoxMode.get())
    +"\n\n#Sets a filter on the KVFinder output, excluding cavities with smaller volumes than this parameter.\n>VOLUME_FILTER\n"+v.get()+"\n\n#Selects surface type to be considered, (0) Solvent Accessible Surface (SAS) or (1) Molecular Surface (VdW).\n>SURF\n"+str(surf.get())
    +"\n\n#This option is used when the user wants to limit the search space around a ligand. (1) for on and (0) for off.\n>LIGAND_ADJUSTMENT\n"+str(ligandAd.get())+"\n\n#If the ligand adjustment is set, this specifies the path for the ligand file.\n>LIGAND_FILE\n"+lf.get(lf.curselection())+"\n\n#Defines the seach radius for the ligand adjustmen mode.\n>DISTANCE_CUTOFF\n"+dc.get())

    f_out.close()


############ Draw box functions ################################

#Deletes the box and all the selections used to buid it

def deleteBox(min_x, max_x, min_y, max_y, min_z, max_z, ang_s, ang2_s):
        global angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6
        angle1 = 0
        angle2 = 0
        lim1 = 0
        lim2 = 0
        lim3 = 0
        lim4 = 0
        lim5 = 0
        lim6 = 0
        cmd.delete("vertices")
        cmd.delete(box_name)

        min_x.set(lim1)
        min_x.configure(state="disabled")
        max_x.set(lim2)
        max_x.configure(state="disabled")
        min_y.set(lim3)
        min_y.configure(state="disabled")
        max_y.set(lim4)
        max_y.configure(state="disabled")
        min_z.set(lim5)
        min_z.configure(state="disabled")
        max_z.set(lim6)
        max_z.configure(state="disabled")
        ang_s.set(angle1)
        ang_s.configure(state="disabled")
        ang2_s.set(angle2)
        ang2_s.configure(state="disabled")


#Each one of the redraw functions is binded with one of the slides bar in the DrawBox Options window. They call the rotateBox function, passing all the parameters used to build a box, with the updated value of one of the parameters

def redrawA1(angle):
        global angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6
        angle1 = (float(angle)/180.0)*3.1415926
        redrawBox(angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6);

def redrawA2(angle):
        global angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6
        angle2 = (float(angle)/180.0)*3.1415926
        redrawBox(angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6);

def redrawL1(lim, min_x, max_x):
        global angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6
        lim1 = float(lim)
        if ((x-lim1) >= (x+lim2)):
            min_x.set(-lim2)
            lim1 = -lim2
        redrawBox(angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6)

def redrawL2(lim, min_x, max_x):
        global angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6
        lim2 = float(lim)
        if ((x-lim1) >= (x+lim2)):
            max_x.set(-lim1)
            lim2 = -lim1
        redrawBox(angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6)

def redrawL3(lim, min_y, max_y):
        global angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6
        lim3 = float(lim)
        if ((y-lim3) >= (y+lim4)):
            min_y.set(-lim4)
            lim3 = -lim4
        redrawBox(angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6)

def redrawL4(lim, min_y, max_y):
        global angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6
        lim4 = float(lim)
        if ((y-lim3) >= (y+lim4)):
            max_y.set(-lim3)
            lim4 = -lim3
        redrawBox(angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6)

def redrawL5(lim, min_z, max_z):
        global angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6
        lim5 = float(lim)
        if ((z-lim5) >= (z+lim6)):
            min_z.set(-lim6)
            lim5 = -lim6
        redrawBox(angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6)

def redrawL6(lim, min_z, max_z):
        global angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6
        lim6 = float(lim)
        if ((z-lim5) >= (z+lim6)):
            max_z.set(-lim5)
            lim6 = -lim5
        redrawBox(angle1, angle2, lim1, lim2, lim3, lim4 ,lim5, lim6)

#Draws a new box, calculating the cordinates around the residues selected.

def drawBox(min_x, max_x, min_y, max_y, min_z, max_z, ang_s, ang2_s, l, exp, par):


        global x, y, z, lim1, lim2, lim3, lim4, lim5, lim6, angle1, angle2

        if not exp.get().replace('.','',1).isdigit():
            tkMessageBox.showerror("Warning", "Invalid padding value", parent=par)
            return


        angle1 = 0
        angle2 = 0
        cmd.delete(box_name)
        selection="(sele)"
        ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)
        x = (minX + maxX)/2
        y = (minY + maxY)/2
        z = (minZ + maxZ)/2
        lim1 = x - (minX - float(exp.get()))
        lim2 = (maxX + float(exp.get())) - x
        lim3 = y - (minY - float(exp.get()))
        lim4 = (maxY + float(exp.get())) - y
        lim5 = z - (minZ - float(exp.get()))
        lim6 = (maxZ + float(exp.get()))  - z


        min_x.configure(state="normal")
        min_x.set(lim1)
        max_x.configure(state="normal")
        max_x.set(lim2)
        min_y.configure(state="normal")
        min_y.set(lim3)
        max_y.configure(state="normal")
        max_y.set(lim4)
        min_z.configure(state="normal")
        min_z.set(lim5)
        max_z.configure(state="normal")
        max_z.set(lim6)
        ang_s.configure(state="normal")
        ang_s.set(angle1)
        ang2_s.configure(state="normal")
        ang2_s.set(angle2)


#Rotates the box. Does the same thing of the drawBox method, but without calculating the initial values of the parameters

def redrawBox(angle1, angle2, lim1, lim2, lim3, lim4, lim5, lim6):



        z1 = lim1*math.sin(float(angle2)) + lim3*math.sin(float(angle1))*math.cos(float(angle2)) - lim5*math.cos(float(angle1))*math.cos(float(angle2)) + z
        y1 = -lim3*math.cos(float(angle1)) + (-lim5)*math.sin(float(angle1)) + y
        x1 = -lim1*math.cos(float(angle2)) - (-lim3)*math.sin(float(angle1))*math.sin(float(angle2)) + (-lim5)*math.cos(float(angle1))*math.sin(float(angle2)) + x

        z2 = -(lim2)*math.sin(float(angle2)) - (-lim3)*math.sin(float(angle1))*math.cos(float(angle2)) + (-lim5)*math.cos(float(angle1))*math.cos(float(angle2)) + z
        y2 = (-lim3)*math.cos(float(angle1)) + (-lim5)*math.sin(float(angle1)) + y
        x2 = (lim2)*math.cos(float(angle2)) - (-lim3)*math.sin(float(angle1))*math.sin(float(angle2)) + (-lim5)*math.cos(float(angle1))*math.sin(float(angle2)) + x

        z3 = -(-lim1)*math.sin(float(angle2)) - (lim4)*math.sin(float(angle1))*math.cos(float(angle2)) + (-lim5)*math.cos(float(angle1))*math.cos(float(angle2)) + z
        y3 = (lim4)*math.cos(float(angle1)) + (-lim5)*math.sin(float(angle1)) + y
        x3 = (-lim1)*math.cos(float(angle2)) - (lim4)*math.sin(float(angle1))*math.sin(float(angle2)) + (-lim5)*math.cos(float(angle1))*math.sin(float(angle2)) + x

        z4 = -(-lim1)*math.sin(float(angle2)) - (-lim3)*math.sin(float(angle1))*math.cos(float(angle2)) + (lim6)*math.cos(float(angle1))*math.cos(float(angle2)) + z
        y4 = (-lim3)*math.cos(float(angle1)) + (lim6)*math.sin(float(angle1)) + y
        x4 = (-lim1)*math.cos(float(angle2)) - (-lim3)*math.sin(float(angle1))*math.sin(float(angle2)) + (lim6)*math.cos(float(angle1))*math.sin(float(angle2)) + x

        z5 = -(lim2)*math.sin(float(angle2)) - (lim4)*math.sin(float(angle1))*math.cos(float(angle2)) + (-lim5)*math.cos(float(angle1))*math.cos(float(angle2)) + z
        y5 = (lim4)*math.cos(float(angle1)) + (-lim5)*math.sin(float(angle1)) + y
        x5 = (lim2)*math.cos(float(angle2)) - (lim4)*math.sin(float(angle1))*math.sin(float(angle2)) + (-lim5)*math.cos(float(angle1))*math.sin(float(angle2)) + x

        z6 = -(lim2)*math.sin(float(angle2)) - (-lim3)*math.sin(float(angle1))*math.cos(float(angle2)) + (lim6)*math.cos(float(angle1))*math.cos(float(angle2)) + z
        y6 = (-lim3)*math.cos(float(angle1)) + (lim6)*math.sin(float(angle1)) + y
        x6 = (lim2)*math.cos(float(angle2)) - (-lim3)*math.sin(float(angle1))*math.sin(float(angle2)) + (lim6)*math.cos(float(angle1))*math.sin(float(angle2)) + x

        z7 = -(-lim1)*math.sin(float(angle2)) - (lim4)*math.sin(float(angle1))*math.cos(float(angle2)) + (lim6)*math.cos(float(angle1))*math.cos(float(angle2)) + z
        y7 = (lim4)*math.cos(float(angle1)) + (lim6)*math.sin(float(angle1)) + y
        x7 = (-lim1)*math.cos(float(angle2)) - (lim4)*math.sin(float(angle1))*math.sin(float(angle2)) + (lim6)*math.cos(float(angle1))*math.sin(float(angle2)) + x

        z8 = -(lim2)*math.sin(float(angle2)) -(lim4)*math.sin(float(angle1))*math.cos(float(angle2)) + (lim6)*math.cos(float(angle1))*math.cos(float(angle2)) + z
        y8 = (lim4)*math.cos(float(angle1)) + (lim6)*math.sin(float(angle1)) + y
        x8 = (lim2)*math.cos(float(angle2)) - (lim4)*math.sin(float(angle1))*math.sin(float(angle2)) + (lim6)*math.cos(float(angle1))*math.sin(float(angle2)) + x



        pymol.stored.list = []
        cmd.iterate(box_name, "stored.list.append((name, color))")
        list_color = pymol.stored.list

        cmd.delete(box_name)

        if len(list_color) > 2:
          for item in list_color:
            at_name = item[0]
            at_c = item[1]
            cmd.set_color(at_name+"color", cmd.get_color_tuple(at_c))
        else:
          for at_name in ["v2","v3","v4","v5","v6","v7","v8","v1x","v1y","v1z","v2x","v3y","v4z"]:
            cmd.set_color(at_name+"color", [0.86,0.86,0.86])



        cmd.pseudoatom(box_name, name="v2", pos=[x2, y2, z2],color="v2color")

        cmd.pseudoatom(box_name, name="v3", pos=[x3, y3, z3],color="v3color")

        cmd.pseudoatom(box_name, name="v4", pos=[x4, y4, z4],color="v4color")

        cmd.pseudoatom(box_name, name="v5", pos=[x5, y5, z5],color="v5color")

        cmd.pseudoatom(box_name, name="v6", pos=[x6, y6, z6],color="v6color")

        cmd.pseudoatom(box_name, name="v7", pos=[x7, y7, z7],color="v7color")

        cmd.pseudoatom(box_name, name="v8", pos=[x8, y8, z8],color="v8color")


        cmd.select("vertices", "(name v3,v7)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v2,v6)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v5,v8)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v2,v5)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v4,v6)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v4,v7)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v3,v5)")
        cmd.bond("vertices", "vertices")

        cmd.select("vertices", "(name v6,v8)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v7,v8)")
        cmd.bond("vertices", "vertices")

        cmd.pseudoatom(box_name, name="v1x", pos=[x1, y1, z1], color = 'red')
        cmd.pseudoatom(box_name, name="v2x", pos=[x2, y2, z2], color = 'red')

        cmd.select("vertices", "(name v1x,v2x)")
        cmd.bond("vertices", "vertices")

        cmd.pseudoatom(box_name, name="v1y", pos=[x1, y1, z1], color = 'forest')
        cmd.pseudoatom(box_name, name="v3y", pos=[x3, y3, z3], color = 'forest')

        cmd.select("vertices", "(name v1y,v3y)")
        cmd.bond("vertices", "vertices")

        cmd.pseudoatom(box_name, name="v4z", pos=[x4, y4, z4], color = 'blue')
        cmd.pseudoatom(box_name, name="v1z", pos=[x1, y1, z1], color = 'blue')
        cmd.select("vertices", "(name v1z,v4z)")
        cmd.bond("vertices", "vertices")
        cmd.delete("vertices")


####################### End of draw box functions ##############################




###################### Listbox refreshers ######################################


#Refreshes the input file listbox in the KVFinder tab

def refreshList(l):
    l.delete(0, l.size())

    for item in cmd.get_names("all"):
        if cmd.get_type(item) == "object:molecule" and item != "box" and item[-16:] != ".KVFinder.output" and item[:7] != 'profile' and item != 'target_exclusive':
            l.insert(l.size(), item)

    l.selection_set(l.size()-1,l.size()-1)


#Refreshes the file selections in the Cavity Profiler tab
#def refreshLists(l_kvp_if, bf):


#    l_kvp_if.delete(0, l_kvp_if.size())
#    bf.delete(0, bf.size())
#    for item in cmd.get_names("all"):
#        if cmd.get_type(item) == "object:molecule" and item != "box" and item[-16:] != ".KVFinder.output" and item[:7] != 'profile' and item != 'target_exclusive':
#            l_kvp_if.insert(l_kvp_if.size(), item)
#            bf.insert(bf.size(), item)

#    l_kvp_if.selection_set(0, l_kvp_if.size())
#    bf.selection_set(0,0)


#Refreshes the results listbox - volume

def ref_kvout_v(o, l_out):
    l_out.delete(0, l_out.size())
    f_aux = open(o.get()+'.KVFinder.results.txt', "r")
    for line in f_aux:
        if line.find(": Volume =") != -1:
            l_out.insert(l_out.size(), line.rstrip('\n').replace(' Volume = ', '').replace(' Angstrons^3', ''))

    f_aux.close()

#Refreshes the results listbox - area

def ref_kvout_a(o, l_out):
    l_out.delete(0, l_out.size())
    f_aux = open(o.get()+'.KVFinder.results.txt', "r")
    for line in f_aux:
        if line.find(": Area =") != -1:
            l_out.insert(l_out.size(), line.rstrip('\n').replace(' Area = ', '').replace(' Angstrons^2', ''))

    f_aux.close()

#Refreshes the residues listbox and the dictionary correlating cavities and residues

def ref_res(o, l_out):
    global residues
    residues = {}
    l_out.delete(0, l_out.size())
    lines = f_aux = open(o.get()+'.KVFinder.results.txt', "r").readlines()
    i = len(lines) - 1
    while (lines[i].find("#Interface") == -1):
        line = lines[i]
        i = i-1
      	if line.split()[0]+" "+line.split()[1] not in residues:
            residues[line.split()[0]+" "+line.split()[1]] = []
       	if len(line.split()) == 4:
	        residues[line.split()[0]+" "+line.split()[1]].append(line.split()[2]+line.split()[3])
        else:
	        residues[line.split()[0]+" "+line.split()[1]].append(line.split()[2])



    for key in sorted(residues.keys()):
        l_out.insert(l_out.size(), key.rstrip(':\n'))



####################### End of the refreshers #########################################

############################ Load Results Function ###################################

def LoadResults(rf, par, l_out, l_out2, l_res,inp_file,lig_file,out_file):
    print rf.get(0,0)
    if (rf.get(0,0).find("KVFinder.results") == -1):
        tkMessageBox.showerror("Warning", "Invalid Results File", parent=par)
        return

    load_file = open(rf.get(0,0),"r")
    for line in load_file:
        if line.find("Input") != -1:
            input_file = line.split()[3]
            if not os.path.isfile(input_file):
                tkMessageBox.showerror("Warning", "Invalid Input File in Results File",parent = par)
                return

        if line.find("Output") != -1:
            output_file = line.split()[3]
            if not os.path.isfile(output_file):
                tkMessageBox.showerror("Warning", "Invalid Output File in Results File",parent = par)
                return


    load_object_input = input_file.split('\\')[len(input_file.split('\\'))-1][:-4]
    load_object_output = output_file.split('\\')[len(output_file.split('\\'))-1][:-4]

    active_objects = cmd.get_names("all")

    i = 0
    j = 1
    k = 1
    while i < len(active_objects):

        if active_objects[i] == load_object_input:
            load_object_input = load_object_input.rstrip(str(j-1))
            load_object_input = load_object_input+str(j)
            j+=1
            i = 0
        if active_objects[i] == load_object_output:
            load_object_output = load_object_output.rstrip(str(k-1))
            load_object_output = load_object_output + str(k)
            k+=1
            i = 0
        i+=1


    cmd.load(input_file, load_object_input)
    cmd.load(output_file, load_object_output)
    out = StringVar()
    out.set(rf.get(0,0)[:-21])

    ref_kvout_v(out, l_out)
    ref_kvout_a(out, l_out2)
    ref_res(out,l_res)

    last_object.set(load_object_output)
    last_input.set(load_object_input)

    inp_file.configure(text="Input File: "+load_object_input)
    lig_file.configure(text='')
    out_file.configure(text="Output File: "+load_object_output)

###################### Results Visual Display Functions ###############################

#selects and highlights cavities selected in the results listbox

def showCavities(event, arg, arg2):

    cavs = []
    for item in arg.curselection():
        cavs.append(arg.get(item)[7:10])
        arg2.selection_set(item,item)

    for item in arg2.curselection():
        if item not in arg.curselection():
            arg2.selection_clear(item,item)

    cmd.set("auto_zoom", 0)
    cmd.delete("cavities")
    cmd.delete("cavities.KVFinder.output")
    if(len(arg.curselection()) < 1):
        return

    control = 0
    for item in cmd.get_names("all"):
        if item == last_object.get():
            control = 1

    if control == 0:
        return

    control = 0
    command = StringVar()
    while len(cavs) > 0:
        command.set(last_object.get()+" and (resname ")
        for item in cavs[:50]:
            command.set(command.get()+item+',')
            cavs.remove(item)

        command.set(command.get()[:-1]+')')
        cmd.select('cavs', command.get())


        if(control == 1):
            cmd.select('cavities','cavs or cavities')
        else:
            cmd.select('cavities', 'cavs')
            control = 1
        cmd.delete('cavs')



    cmd.create("cavities.KVFinder.output", "cavities")
    cmd.delete("cavities")
    if(ChargeDisplay.get() == 0):
        cmd.color("blue", "cavities.KVFinder.output")
    else:
        cmd.do('spectrum pc, red_white_blue, cavities.KVFinder.output, minimum=-3, maximum=3')
    cmd.disable(last_object.get())
    cmd.enable(last_object.get())

    cmd.set("auto_zoom", 1)



#selects and highlights cavities selected in the results listbox

def showCavities2(event, arg, arg2):

    cavs = []
    for item in arg.curselection():
        cavs.append(arg.get(item)[7:10])
        arg2.selection_set(item,item)

    for item in arg2.curselection():
        if item not in arg.curselection():
            arg2.selection_clear(item,item)

    cmd.set("auto_zoom", 0)
    cmd.delete("cavities")
    cmd.delete("cavities.KVFinder.output")
    if(len(arg.curselection()) < 1):
        return

    control = 0
    for item in cmd.get_names("all"):
        if item == last_object.get():
            control = 1

    if control == 0:
        return

    control = 0
    command = StringVar()
    while len(cavs) > 0:
        command.set(last_object.get()+" and (resname ")
        for item in cavs[:50]:
            command.set(command.get()+item+',')
            cavs.remove(item)

        command.set(command.get()[:-1]+')')
        cmd.select('cavs', command.get())


        if(control == 1):
            cmd.select('cavities','cavs or cavities')
        else:
            cmd.select('cavities', 'cavs')
            control = 1
        cmd.delete('cavs')



    cmd.create("cavities.KVFinder.output", "cavities")
    cmd.delete("cavities")
    if(ChargeDisplay.get() == 0):
        cmd.color("blue", "cavities.KVFinder.output")
    else:
        cmd.do('spectrum pc, red_white_blue, cavities.KVFinder.output, minimum=-3, maximum=3')
    cmd.disable(last_object.get())
    cmd.enable(last_object.get())

    cmd.set("auto_zoom", 1)

#selects and highlights residues of the cavities selected in the results listbox

def showResidues(event, arg):

    chain = ''

    res_lis = {}

    for item in arg.curselection():
        for x in residues[arg.get(item)]:
            if x[-1:] not in res_lis:
                res_lis[x[-1:]] = []
                res_lis[x[-1:]].append(x[:-1])
            else:
                res_lis[x[-1:]].append(x[:-1])

    cmd.set("auto_zoom", 0)
    cmd.delete("residues")
    cmd.delete("residues.KVFinder.output")
    if(len(arg.curselection()) < 1):
        return

    control = 0
    for item in cmd.get_names("all"):
        if item == last_object.get():
            control = 1

    if control == 0:
        return

    control = 0
    cont = 0
    command = StringVar()
    for item in res_lis:
            while len(res_lis[item]) > 0:
                command.set(last_input.get()+" and chain "+item+" and (resid ")
                for x in res_lis[item][:150]:
                    command.set(command.get()+x+',')
                    res_lis[item].remove(x)



                command.set(command.get()[:-1]+')')
                cmd.select('res', command.get())



                if(control == 1):
                    cmd.select('residues','res or residues')
                else:
                    cmd.select('residues', 'res')
                    control = 1

                cmd.delete('res')



    cmd.create("residues.KVFinder.output", "residues")
    cmd.delete("residues")
    cmd.show('sticks', 'residues.KVFinder.output')
    cmd.disable(last_object.get())
    cmd.enable(last_object.get())

    cmd.set("auto_zoom", 1)



##################### End of the Visual Display ############################

########## Functions to execute KVFinder and Cavity Profiler ###############


#Function that saves the Parameters.txt and runs the KVFinder
def runKVFinder(d, s, p, v, surf, lf, dc, l, o, l_out, p_out,inp_file,lig_file,out_file,l_res,l_out2, par,par2):

    if( lim1==0 and lim2==0 and lim3==0 and lim4==0 and lim5==0 and lim6==0 and WholeProt.get()==0):
        tkMessageBox.showerror("Warning", "Search space not defined. Please, define a search space by selecting some residues and clicking in Draw New Box or use the Whole Protein mode", parent=par)
        return

    if not s.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid step size", parent=par)
        return

    if not p.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid probe size", parent=par)
        return

    if not p_out.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid probe out size", parent=par)
        return

    if (WholeProt.get() == 1 or BoxMode.get() == 1) and StepRed.get() == 0 and float(s.get()) <= 1 and  float(p_out.get()) >= 20:
        if not tkMessageBox.askyesno("Warning", "Given the Probe Out defined, this operation may take a really long time. Should it proceed?", parent=par):
            return

    if not v.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid minimum cavity volume", parent=par)
        return

    if not dc.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid distance cutoff", parent=par)
        return

    if not os.path.isfile(d.get()):
        tkMessageBox.showerror("Warning", "Dictionary Filename", parent=par)
        return

    if not l.curselection():
        tkMessageBox.showerror("Warning", "Select a input file", parent=par)
        return

    if len(o.get()) == 0:
         tkMessageBox.showerror("Warning", "Output file not defined", parent=par)
         return

    if ligandAd.get() == 1 and not lf.curselection():
        tkMessageBox.showerror("Warning", "Ligand Adjustment selected with invalid ligand filename", parent=par)
        return

    control = 1
    if ligandAd.get() == 1:
        for x in cmd.get_names("all"):
            if x == lf.get(lf.curselection()):
                control = 0

        if os.path.isfile(lf.get(lf.curselection())) == 1:
            control = 0

        if control == 1:
            tkMessageBox.showerror("Warning", "Ligand Adjustment selected with invalid ligand", parent=par)
            return

    if ligandAd.get() == 1 and float(dc.get()) <= 0:
        tkMessageBox.showerror("Warning", "Ligand Adjustment selected with invalid Distance Cutoff", parent=par)
        return

    if WholeProt.get() == 1 and (float(p.get()) >= float(p_out.get())):
        tkMessageBox.showerror("Warning", "Invalid probe sizes for Whole protein mode. Probe Out must be bigger than probe", parent=par)
        return

    if BoxMode.get() == 1 and (float(p.get()) >= float(p_out.get())):
        tkMessageBox.showerror("Warning", "Invalid probe sizes for Probe Out Adjustment mode. Probe Out must be bigger than probe", parent=par)
        return

    control = 1
    for x in cmd.get_names("all"):
        if x == l.get(l.curselection()):
            control = 0

    if control == 1:
        tkMessageBox.showerror("Warning", "Invalid input file", parent=par)
        return

    cmd.delete("cavities.KVFinder.output")
    cmd.delete("residues.KVFinder.output")

    tmpLigFile = StringVar()
    tmpLigFile.set(lf.get(lf.curselection()))

    if os.path.isdir(path.get()+"KV_Files") == 0:
        os.makedirs(path.get()+'KV_Files')

    for x in cmd.get_names("all"):
        if lf.get(lf.curselection()) == x:
            cmd.save(path.get()+'KV_Files\KVFinder_input_'+lf.get(lf.curselection())+'.pdb', lf.get(lf.curselection()), 0, "pdb")
            tmpLigFile.set(path.get()+'KV_Files\KVFinder_input_'+lf.get(lf.curselection())+'.pdb')


    i = len(o.get()) - 1
    while o.get()[i] != '\\' and i > -1:
        i = i - 1
    i = i + 1

    for x in cmd.get_names("all"):
        if o.get()[i:]+".KVFinder.output" == x:
            if(OverW.get() == 0):
                tkMessageBox.showerror("Warning", "This operation will create a object with the same name of an existing object.\nChange the output file or delete the existing object.", parent=par)
                return
            else:
                cmd.delete(x)



    cmd.save(path.get()+'KV_Files\\'+o.get()+'.KVFinder.input.pdb', l.get(l.curselection()), 0, "pdb")
    infile.set(path.get()+'KV_Files\\'+o.get()+'.KVFinder.input.pdb')

    out_complete = StringVar()
    out_complete.set(path.get()+o.get())

    saveBox(d, s, p, v, surf, tmpLigFile, dc, out_complete,p_out,par2)


    results = open(o.get()+".KVFinder.results.txt","w+")
    results.write("KVFinder Input = "+'KV_Files\\'+o.get()+'.KVFinder.input.pdb\nKVFinder Output = '+o.get()+".KVFinder.output.pdb\n\n#KVFinder Results:\n")
    results.close()


    if os.path.isfile(path.get()+'KV_Files\\kvfinder_error.log') == 1:
        os.system("del "+'KV_Files\\kvfinder_error.log')

    print "Running KVFinder for "+infile.get()+"..."
    os.system("\""+os.getenv('KVFinder_PATH')+'\\KVFinder.exe\" 2> '+'KV_Files\\kvfinder_error.log 1> '+o.get()+'.KVFinder.output.tmp')
    print "done!"

    log = open(path.get()+"KV_Files\\KVFinder.log","a+")
    tmp = open(path.get()+o.get()+'.KVFinder.output.tmp',"r")

    for line in tmp:
        log.write(line)
    tmp.close()
    log.close()

    tmp = open(path.get()+o.get()+'.KVFinder.output.tmp',"r")
    results = open(o.get()+".KVFinder.results.txt","a+")

    for line in tmp:
        if line.find("Cavity")!=-1:
            results.write(line)

    results.close()
    tmp.close()



    os.system('del '+o.get()+'.KVFinder.output.tmp')
    last_object.set(o.get()[i:]+'.KVFinder.output')

    f_in = open(path.get()+'KV_Files\kvfinder_error.log', 'r')
    for line in f_in:
        if line[:18] == "Segmentation fault":
            tkMessageBox.showerror("Warning", "There was an error in the KVFinder software. Try using a smaller box or a bigger step size.", parent=par)
            os.system("del "+'KV_Files\kvfinder_error.log')
            f_in.close()
            return
    f_in.close()
    os.system("del "+'KV_Files\kvfinder_error.log')



    inp_file.configure(text="Input File: "+l.get(l.curselection()))
    if ligandAd.get() == 1:
        if len(lf.get(lf.curselection())) > 30:
            j = 0
            while lf.get(lf.curselection())[-30:][j] != '/' and j < len(lf.get(lf.curselection())[-30:]) - 1:
                j = j + 1
            if len(lf.get(lf.curselection())[-30:]) - j == 1:
                j = 0
            lig_file.configure(text="Ligand File: ..."+lf.get(lf.curselection())[-(30-j):])
        else:
            lig_file.configure(text="Ligand File: "+lf.get(lf.curselection()))
    else:
        lig_file.configure(text='')
    out_file.configure(text="Output File: "+last_object.get())


    if os.path.isfile(path.get()+o.get()+".KVFinder.output.pdb") == 1:
        cmd.set("auto_zoom",0)
        cmd.load(path.get()+o.get()+".KVFinder.output.pdb")
        cmd.set("auto_zoom",1)

    os.system('del Parameters.txt')
    os.system('copy KVParameters.txt '+path.get()+'KV_Files\\KVParameters_'+o.get()[i:])

    results = open(o.get()+".KVFinder.results.txt","a")
    results.write("\n#Interface Residues for Each Cavity\n")

    cavres = open("cavres.txt", "r")
    for line in cavres:
        results.write(line)

    results.close()
    cavres.close()
    os.system("del cavres.txt")

    ref_kvout_v(out_complete, l_out)
    ref_kvout_a(out_complete, l_out2)
    ref_res(out_complete,l_res)

    last_input.set(l.get(l.curselection()))


#Function that saves the Parameters.txt and runs the KVFinder
def runKVFinderBasic(d, s, p, v, o, surf, lf, dc, l, l_out, p_out,inp_file,lig_file,out_file,l_res,l_out2, par,par2):

    if not s.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid step size", parent=par)
        return

    if not p.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid probe size", parent=par)
        return

    if not p_out.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid probe out size", parent=par)
        return

    if float(s.get()) <= 1 and  float(p_out.get()) >= 30:
        if not tkMessageBox.askyesno("Warning", "Given the Probe Out defined, this operation may take a really long time. Should it proceed?", parent=par):
            return


    if not os.path.isfile(d.get()):
        tkMessageBox.showerror("Warning", "Invalid Dictionary Filename. Check the configuration tab.", parent=par)
        return


    if not v.get().replace('.','',1).isdigit():
        tkMessageBox.showerror("Warning", "Invalid minimum cavity volume", parent=par)
        return

    if not l.curselection():
        tkMessageBox.showerror("Warning", "Select a input file", parent=par)
        return


    if (float(p.get()) >= float(p_out.get())):
        tkMessageBox.showerror("Warning", "Invalid probe sizes for Whole protein mode. Probe Out must be bigger than probe", parent=par)
        return

    control = 1
    for x in cmd.get_names("all"):
        if x == l.get(l.curselection()):
            control = 0

    if control == 1:
        tkMessageBox.showerror("Warning", "Invalid input file", parent=par)
        return

    cmd.delete("cavities.KVFinder.output")
    cmd.delete("residues.KVFinder.output")

    if os.path.isdir(path.get()+"KV_Files") == 0:
        os.makedirs(path.get()+'KV_Files')


    cmd.save(path.get()+'KV_Files\\'+o.get()+'.KVFinder.input.pdb', l.get(l.curselection()), 0, "pdb")
    infile.set(path.get()+'KV_Files\\'+o.get()+'.KVFinder.input.pdb')

    surfTmp = BooleanVar()
    surfTmp.set(1)

    StepRedTmp = BooleanVar()
    StepRedTmp.set(StepRed.get())
    StepRed.set(0)

    BoxModeTmp = BooleanVar()
    BoxModeTmp.set(BoxMode.get())
    BoxMode.set(0)

    WholeProtTmp = BooleanVar()
    WholeProtTmp.set(WholeProt.get())
    WholeProt.set(1)

    ligandAdTmp = BooleanVar()
    ligandAdTmp.set(ligandAd.get())
    ligandAd.set(0)

    lf_tmp = StringVar()
    lf_tmp.set('-')

    saveBox(d, s, p, v, surfTmp, lf_tmp, dc, o ,p_out,par2)

    WholeProt.set(WholeProtTmp.get())
    BoxMode.set(BoxModeTmp.get())
    StepRed.set(StepRedTmp.get())
    ligandAd.set(ligandAdTmp.get())

    results = open(o.get()+".KVFinder.results.txt","w+")
    results.write("KVFinder Input = "+'KV_Files\\'+o.get()+'.KVFinder.input.pdb\nKVFinder Output = '+o.get()+".KVFinder.output.pdb\n\n#KVFinder Results:\n")
    results.close()

    if os.path.isfile(path.get()+'KV_Files\\kvfinder_error.log') == 1:
        os.system("del "+'KV_Files\\kvfinder_error.log')
    print "Running KVFinder for "+infile.get()+"..."
    os.system("\""+os.getenv('KVFinder_PATH')+'\\KVFinder.exe\" 2> '+'KV_Files\\kvfinder_error.log 1> '+o.get()+'.KVFinder.output.tmp')
    print "done!"

    log = open("KV_Files\\KVFinder.log","a+")
    tmp = open(o.get()+'.KVFinder.output.tmp',"r")

    for line in tmp:
        log.write(line)
    tmp.close()
    log.close()

    tmp = open(o.get()+'.KVFinder.output.tmp',"r")
    results = open(o.get()+".KVFinder.results.txt","a+")

    for line in tmp:
        if line.find("Cavity")!=-1:
            results.write(line)

    results.close()
    tmp.close()



    os.system('del '+o.get()+'.KVFinder.output.tmp')
    last_object.set(o.get()+'.KVFinder.output')

    f_in = open(path.get()+'KV_Files\kvfinder_error.log', 'r')
    for line in f_in:
        if line[:18] == "Segmentation fault":
            tkMessageBox.showerror("Warning", "There was an error in the KVFinder software. Try using a smaller box or a bigger step size.", parent=par)
            os.system("del "+'KV_Files\kvfinder_error.log')
            f_in.close()
            return
    f_in.close()
    os.system("del "+'KV_Files\kvfinder_error.log')



    inp_file.configure(text="Input File: "+l.get(l.curselection()))
    last_input.set(l.get(l.curselection()))
    if ligandAd.get() == 1:
        if len(lf.get(lf.curselection())) > 30:
            j = 0
            while lf.get(lf.curselection())[-30:][j] != '\\' and j < len(lf.get(lf.curselection())[-30:]) - 1:
                j = j + 1
            if len(lf.get(lf.curselection())[-30:]) - j == 1:
                j = 0
            lig_file.configure(text="Ligand File: ..."+lf.get(lf.curselection())[-(30-j):])
        else:
            lig_file.configure(text="Ligand File: "+lf.get(lf.curselection()))
    else:
        lig_file.configure(text='')
    out_file.configure(text="Output File: "+last_object.get())


    if os.path.isfile(o.get()+".KVFinder.output.pdb") == 1:
        cmd.set("auto_zoom",0)
        cmd.load(o.get()+".KVFinder.output.pdb")
        cmd.set("auto_zoom",1)

    os.system('del Parameters.txt')
    os.system('copy KVParameters.txt KV_File\\KVParameters_'+o.get())

    results = open(o.get()+".KVFinder.results.txt","a")
    results.write("\n#Interface Residues for Each Cavity\n")

    cavres = open("cavres.txt", "r")
    for line in cavres:
        results.write(line)

    results.close()
    cavres.close()
    os.system("del cavres.txt")

    ref_kvout_v(o, l_out)
    ref_kvout_a(o, l_out2)
    ref_res(o,l_res)




############### End of the exection functions ######################



############# Auxiliar Functions ########################


#Opens a selection file window, used to select a ligand file in the Ligand Adjustment option

def selectFile(lf):

    if lf.get(0) == '-':
        lf.delete(0,0)

    filename = tkFileDialog.askopenfilename()
    lf.insert(lf.size(),filename)
    lf.selection_clear(0, lf.size()-1)
    lf.selection_set(lf.size()-1,lf.size()-1)


#Opens a selection file window, used to select a ligand file in the Load Results Option

def selectFileResults(rf):

    rf.delete(0,0)

    filename = tkFileDialog.askopenfilename(filetypes=[("KVFinder Results File","*KVFinder.results.txt")])
    rf.insert(rf.size(),filename)
    rf.selection_clear(0, rf.size()-1)
    rf.selection_set(rf.size()-1,rf.size()-1)

#Opens a selection folder window, used to select the output folder

def selectFolder(folder):
    oldFolder = folder.get()

    folder.configure(state='normal')
    folder.delete(0,len(folder.get()))

    folderName = tkFileDialog.askdirectory(mustexist=1)
    if len(folderName) == 0:
        folderName = oldFolder

    folder.insert(0,folderName)
    folder.configure(state='readonly')

    path.set(folder.get()+'\\')


#restore congifuration tab

def restoreConfig(d):

    d.delete(0, len(d.get()))
    d.insert(0, os.getenv('KVFinder_PATH')+'\\dictionary')
    d.xview(d.index('end')-20)


def restoreOptions(s, p, v, sas, ms, checkLA, checkSR, lf, dc, o, p_out, checkOW, checkWP, checkBM, exp, surf):

    global StepRed, ligandAd

    s.delete(0, len(s.get()))
    s.insert(0, "0.6")

    p.delete(0, len(p.get()))
    p.insert(0, "1.4")

    exp.delete(0, len(exp.get()))
    exp.insert(0, "3.5")

    p_out.delete(0, len(p_out.get()))
    p_out.insert(0, "4.0")

    v.delete(0, len(v.get()))
    v.insert(0, "5.0")


    sas.deselect()
    ms.invoke()
    surf.set(1)

    checkLA.deselect()
    checkSR.select()

    StepRed.set(1)
    ligandAd.set(0)

    dc.delete(0, len(dc.get()))
    dc.insert(0, "0.0")

    lf.delete(0, lf.size())
    lf.insert(lf.size(), '-')
    lf.selection_set(0,0)

    o.xview(len(o.get())-1)


    WholeProt.set(0)
    OverW.set(0)
    ChargeDisplay.set(0)
    BoxMode.set(1)


    checkOW.deselect()
    checkWP.deselect()
    checkBM.select()

#restores the default values for most of the KVFinder options

def restore(d, s, p, v, sas, ms, checkLA, checkSR, lf, dc, o, p_out, checkOW, checkWP, checkBM, exp, s_basic, p_basic, p_out_basic, surf, v_basic):

    global StepRed, ligandAd, OverW, chargeDisplay, WholeProt, BoxMode

    d.delete(0, len(d.get()))
    d.insert(0, os.getenv('KVFinder_PATH')+'\\dictionary')
    d.xview(d.index('end')-20)


    s.delete(0, len(s.get()))
    s.insert(0, "0.6")

    s_basic.delete(0, len(s_basic.get()))
    s_basic.insert(0, "0.6")

    p.delete(0, len(p.get()))
    p.insert(0, "1.4")
    p_basic.delete(0, len(p_basic.get()))
    p_basic.insert(0, "1.4")

    exp.delete(0, len(exp.get()))
    exp.insert(0, "3.5")

    p_out.delete(0, len(p_out.get()))
    p_out.insert(0, "4.0")
    p_out_basic.delete(0, len(p_out_basic.get()))
    p_out_basic.insert(0, "4.0")


    v.delete(0, len(v.get()))
    v.insert(0, "5.0")

    v_basic.delete(0, len(v_basic.get()))
    v_basic.insert(0, "5.0")

    sas.deselect()
    ms.invoke()
    surf.set(1)

    checkLA.deselect()
    checkSR.select()

    StepRed.set(1)
    ligandAd.set(0)

    dc.delete(0, len(dc.get()))
    dc.insert(0, "0.0")

    lf.delete(0, lf.size())
    lf.insert(lf.size(), '-')
    lf.selection_set(0,0)

    o.xview(len(o.get())-1)

    WholeProt.set(0)
    OverW.set(0)

    ChargeDisplay.set(0)
    BoxMode.set(1)

    checkOW.deselect()
    checkWP.deselect()
    checkBM.select()


############## Interface code ########################

# tk interpretator for the plugin
# Creates the main window
def testDialog():

##Check if the system variable is defined. It won't work in super user mode
    global outfile
    global path
    refList = lambda: refreshList(l)
    refListBasic = lambda: refreshList(l_basic)


    if(os.getenv("KVFinder_PATH") ==None):
        tkMessageBox.showerror("Warning", "The system variable KVFinder_PATH is not defined. Please read the KVFinder Installation Guide for more information.")
        return

    root = Tk() # start it

    root.title('KVFinder') # title
    mainframe  = Frame(root)
    mainframe.grid()

    advanced_frame = Frame(mainframe)

    # if resize window, keep all in place.
    mainframe.columnconfigure(0, weight=1)
    mainframe.rowconfigure(0, weight=1)

    notebook = Pmw.NoteBook(mainframe)
    notebook.grid()

    basic = notebook.add('Basic')
    notebook.tab('Basic').focus_set()

    advanced_tab = notebook.add('Advanced')


    notebook_adv = Pmw.NoteBook(advanced_tab)
    notebook_adv.grid()

    main = notebook.add('Results')

    config = notebook.add('Configurations')

    about = notebook.add('About')

    drbox = notebook_adv.add('Main')
    opt = notebook_adv.add('Options')

################### Basic ###########################

    Label(basic, text = "").pack(pady=3)

    l_frame_basic = Frame(basic)
    scrollbar_b = Scrollbar(l_frame_basic, orient="vertical")
    Label(l_frame_basic, text = "Input File:").pack()
    l_basic = Listbox(l_frame_basic, height = 12, width=30, yscrollcommand=scrollbar_b.set,exportselection=0)

    scrollbar_b.config(command=l_basic.yview)
    scrollbar_b.pack(side=RIGHT, fill=Y)

    l_basic.pack(side=LEFT,fill=BOTH, expand=1)

    l_frame_basic.pack(pady=2)


    Button(basic, text="Refresh List", command = refListBasic).pack()
    Label(basic, text = "").pack(pady=4)


    kvfinderSetBasic = lambda: runKVFinderBasic(d, s_basic, p_basic, v_basic, baseName_basic, surf, lf, dc, l_basic, l_out, p_out_basic,inp_file,lig_file,out_file,l_res,l_out2, basic,right_drbox)
    Button(basic, text="Run KVFinder", fg='DarkGreen', command = kvfinderSetBasic).pack()


    Label(basic, text = "").pack(pady=12)
    Label(basic, text = "Basic Options:").pack(pady=3)
    Label(basic, text = "").pack(pady=6)

    Label(basic, text = "Step Size ("+unichr(0x212b)+"):").pack(pady=3)
    s_basic = Entry(basic,  takefocus=0, width=30)
    s_basic.pack(pady=3)

    Label(basic, text = "Probe In Size ("+unichr(0x212b)+"):").pack(pady=3)
    p_basic = Entry(basic,  takefocus=0, width=30)
    p_basic.pack(pady=3)

    Label(basic, text = "Probe Out Size ("+unichr(0x212b)+"):").pack(pady=3)
    p_out_basic = Entry(basic, takefocus=0, width=30)
    p_out_basic.pack(pady=3)

    Label(basic, text = "Minimum Cavity Volume ("+unichr(0x212b)+unichr(0x00b3)+"):").pack(pady=3)
    v_basic = Entry(basic, takefocus=0, width=30)
    v_basic.pack(pady=3)

    Label(basic, text = "Output Base Name:").pack(pady=3)
    baseName_basic = Entry(basic, takefocus=0, width=30)
    baseName_basic.insert(0,outfileBasic.get())
    baseName_basic.xview(baseName_basic.index('end')-20)
    baseName_basic.pack(pady=3)

    checkOW = Checkbutton(basic, text="Overwrite Existing File", command=c_OW, takefocus=0)
    checkOW.pack()


################### End Basic #####################

#################### Results #######################

    left_main = Frame(main)
    right_main = Frame(main)

    left_main.pack(side='left', fill='both', padx = 25)
    right_main.pack(side='left', fill='both', padx = 25)

    #Label(main, text="").pack(pady=3)
    info_frame = Frame(left_main)
    inp_file = Label(info_frame, text="")
    inp_file.grid(column=1,row=1, stick='W')
    lig_file = Label(info_frame, text="")
    lig_file.grid(column=1,row=2,stick='W')
    out_file = Label(info_frame, text="")
    out_file.grid(column=1,row=3,stick='W')
    info_frame.pack()
    Label(left_main, text="").pack(pady=3)
    out_frame = Frame(left_main)
    scrollbar_of = Scrollbar(out_frame, orient="vertical")
    scrollbar_ofh = Scrollbar(out_frame, orient="horizontal")
    Label(out_frame, text = "KVFinder Output - Volume ("+unichr(0x212b)+unichr(0x00b3)+"):").pack()
    l_out = Listbox(out_frame, height = 10, width=30, yscrollcommand=scrollbar_of.set, xscrollcommand=scrollbar_ofh.set,selectmode='multiple', exportselection=0)

    scrollbar_of.config(command=l_out.yview)
    scrollbar_of.pack(side=RIGHT, fill=Y)
    scrollbar_ofh.config(command=l_out.xview)
    scrollbar_ofh.pack(side=BOTTOM, fill=X)
    l_out.pack(side=LEFT,fill=BOTH, expand=1)


    out_frame_a = Frame(left_main)
    scrollbar_of_a = Scrollbar(out_frame_a, orient="vertical")
    scrollbar_ofh_a = Scrollbar(out_frame_a, orient="horizontal")
    Label(out_frame_a, text = "KVFinder Output - Area("+unichr(0x212b)+unichr(0x00b2)+"):").pack()
    l_out2 = Listbox(out_frame_a, height = 10, width=30, yscrollcommand=scrollbar_of_a.set, xscrollcommand=scrollbar_ofh_a.set,selectmode='multiple', exportselection=0)

    scrollbar_of_a.config(command=l_out2.yview)
    scrollbar_of_a.pack(side=RIGHT, fill=Y)
    scrollbar_ofh_a.config(command=l_out2.xview)
    scrollbar_ofh_a.pack(side=BOTTOM, fill=X)
    l_out2.pack(side=LEFT,fill=BOTH, expand=1)


    l_out.bind('<ButtonRelease-1>', lambda event, arg=l_out, arg2=l_out2: showCavities(event, arg, arg2))
    l_out2.bind('<ButtonRelease-1>', lambda event, arg=l_out2, arg2=l_out: showCavities2(event, arg, arg2))


    out_frame.pack()
    Label(left_main, text="").pack(pady=2)
    out_frame_a.pack()


#    Label(left_main, text="").pack(pady=2)
#    checkCD = Checkbutton(left_main, text="Color by charge", command=c_CD, takefocus=0)
#    checkCD.pack()

    Label(right_main, text="").pack(pady=30)
    #Label(right_main, text="").pack(pady=2)
    res_frame = Frame(right_main)
    scrollbar_res = Scrollbar(res_frame, orient="vertical")
    Label(res_frame, text = "Residues Around:").pack()
    l_res = Listbox(res_frame, height = 25, width=30, yscrollcommand=scrollbar_res.set, selectmode='multiple', exportselection=0)

    scrollbar_res.config(command=l_res.yview)
    scrollbar_res.pack(side=RIGHT, fill=Y)
    l_res.pack(side=LEFT,fill=BOTH, expand=1)


    l_res.bind('<ButtonRelease-1>', lambda event, arg=l_res: showResidues(event, arg))

   #Load Results

    Label(left_main, text = "").pack(pady="3")

    Label(left_main, text = "Load Results:").pack()


    rf = Listbox(left_main, height = 1, width=30, exportselection=0)

    #sb_opt_lf.config(command=lf.yview)
    #sb_opt_lf.pack(side=RIGHT, fill=Y)

    rf.pack()

    selFileResults = lambda: selectFileResults(rf)
    setLoadResults = lambda: LoadResults(rf, left_main, l_out, l_out2, l_res, inp_file, lig_file, out_file)
    results_buttons = Frame(left_main)

    Label(left_main, text = "").pack(pady="3")

    Button(results_buttons, text = 'Select File', command = selFileResults).pack(side=LEFT)
    Button(results_buttons,  text = 'Load Results', command = setLoadResults).pack(side=RIGHT)
    results_buttons.pack()
    res_frame.pack()



#########################End of Main ########################




################ Constructs the Drawbox Frame ##########################
    left_drbox = Frame(drbox)
    right_drbox = Frame(drbox)

    left_drbox.pack(side='left', fill='both', padx = 15)
    right_drbox.pack(side='left', fill='both', padx = 15)

    Label(left_drbox, text = "").pack(pady=3)
    Label(left_drbox, text = 'Select a few residues and press Draw New Box').pack()

    drwBut = Frame(left_drbox)
    drawBoxTab = lambda: drawBox(min_x, max_x, min_y, max_y, min_z, max_z, ang_s, ang2_s, l, exp, left_drbox)
    DNBButton = Button(drwBut, text="Draw New Box", command = drawBoxTab)
    DNBButton.pack(side='left')
    deleteBoxTab = lambda: deleteBox(min_x, max_x, min_y, max_y, min_z, max_z, ang_s, ang2_s)
    Button(drwBut, text="Delete Box", command = deleteBoxTab).pack(side = 'left')
    drwBut.pack(pady = 3)

    Label(left_drbox, text = "").pack(pady=1)


    Label(left_drbox, text = "Box Padding:").pack()
    exp = Entry(left_drbox, takefocus=0, width=30)
    exp.pack(pady=2)


    Label(left_drbox, text = "").pack(pady=1)

    redrawL1Set = lambda(lim): redrawL1(lim, min_x, max_x)

    limA = DoubleVar()
    limA.set(lim1)


    minxLabel = Frame(left_drbox)

    Label(minxLabel, text='Minimum').grid(column=0, row=2)
    Label(minxLabel, text='X', fg='red').grid(column=1, row=2)
    Label(minxLabel, text=':  ').grid(column=2, row=2)


    minxLabel.pack()
    min_x = Scale(minxLabel, resolution = 0.1, orient=HORIZONTAL, length = 150, from_=-50, to = 50, variable = limA, command = redrawL1Set, state="disabled")
    min_x.grid(column=3, row=0, rowspan=3)

    redrawL2Set = lambda(lim): redrawL2(lim, min_x, max_x)

    limB = DoubleVar()
    limB.set(lim2)


    maxxLabel = Frame(left_drbox)

    Label(maxxLabel, text='Maximum').grid(column=0, row=2)
    Label(maxxLabel, text='X', fg='red').grid(column=1, row=2)
    Label(maxxLabel, text=': ').grid(column=2, row=2)


    maxxLabel.pack()




    max_x = Scale(maxxLabel, resolution = 0.1, orient=HORIZONTAL, length = 150, from_=-50, to = 50, variable = limB, command = redrawL2Set, state="disabled")
    max_x.grid(column=3, row=0, rowspan=3)

    redrawL3Set = lambda(lim): redrawL3(lim, min_y, max_y)

    limC = DoubleVar()
    limC.set(lim3)

    minyLabel = Frame(left_drbox)

    Label(minyLabel, text='Minimum').grid(column=0, row=2)
    Label(minyLabel, text='Y', fg='dark green').grid(column=1, row=2)
    Label(minyLabel, text=':  ').grid(column=2, row=2)


    minyLabel.pack()


    min_y = Scale(minyLabel, resolution = 0.1, orient=HORIZONTAL, length = 150, from_=-50, to = 50, variable = limC, command = redrawL3Set, state="disabled")
    min_y.grid(column=3, row=0, rowspan=3)

    redrawL4Set = lambda(lim): redrawL4(lim, min_y, max_y)

    limD = DoubleVar()
    limD.set(lim4)

    maxyLabel = Frame(left_drbox)

    Label(maxyLabel, text='Maximum').grid(column=0, row=2)
    Label(maxyLabel, text='Y', fg='dark green').grid(column=1, row=2)
    Label(maxyLabel, text=': ').grid(column=2, row=2)


    maxyLabel.pack()


    max_y = Scale(maxyLabel, resolution = 0.1, orient=HORIZONTAL, length = 150, from_=-50, to = 50, variable = limD, command = redrawL4Set, state="disabled")
    max_y.grid(column=3, row=0, rowspan=3)

    redrawL5Set = lambda(lim): redrawL5(lim, min_z, max_z)

    limE = DoubleVar()
    limE.set(lim5)

    minzLabel = Frame(left_drbox)

    Label(minzLabel, text='Minimum').grid(column=0, row=2)
    Label(minzLabel, text='Z', fg='blue').grid(column=1, row=2)
    Label(minzLabel, text=':  ').grid(column=2, row=2)


    minzLabel.pack()


    min_z = Scale(minzLabel, resolution = 0.1, orient=HORIZONTAL, length = 150, from_=-50, to = 50, variable = limE, command = redrawL5Set, state="disabled")
    min_z.grid(column=3, row=0, rowspan=3)

    redrawL6Set = lambda(lim): redrawL6(lim, min_z, max_z)

    limF = DoubleVar()
    limF.set(lim6)

    maxzLabel = Frame(left_drbox)

    Label(maxzLabel, text='Maximum').grid(column=0, row=2)
    Label(maxzLabel, text='Z', fg='blue').grid(column=1, row=2)
    Label(maxzLabel, text=': ').grid(column=2, row=2)


    maxzLabel.pack()

    max_z = Scale(maxzLabel, resolution = 0.1, orient=HORIZONTAL, length = 150, from_=-50, to = 50, variable = limF, command = redrawL6Set, state="disabled")
    max_z.grid(column=3, row=0, rowspan=3)



    ang = DoubleVar()
    ang.set(float(angle1)*180/3.1415926)

    ang1Label = Frame(left_drbox)
    Label(ang1Label, text='Angle:           ').grid(column=0, row=2)
    ang_s = Scale(ang1Label, resolution = 0.1, orient=HORIZONTAL, length = 150, from_=-180, to = 180, variable = ang, command = redrawA1, state = "disabled")
    ang_s.grid(column=3, row=0, rowspan=3)
    ang1Label.pack()

    ang2 = DoubleVar()
    ang2.set(float(angle2)*180/3.1415926)

    ang2Label = Frame(left_drbox)
    Label(ang2Label, text='Angle2:          ').grid(column=0, row=2)
    ang2_s = Scale(ang2Label, resolution = 0.1, orient=HORIZONTAL, length = 150, from_=-180, to = 180, variable = ang2, command = redrawA2, state = "disabled")
    ang2_s.grid(column=3, row=0, rowspan=3)
    ang2Label.pack()


    Label(right_drbox, text = "").pack(pady=15)
    l_frame = Frame(right_drbox)
    scrollbar = Scrollbar(l_frame, orient="vertical")
    Label(l_frame, text = "Input File:").pack()
    l = Listbox(l_frame, height = 11, width=30, yscrollcommand=scrollbar.set,exportselection=0)

    scrollbar.config(command=l.yview)
    scrollbar.pack(side=RIGHT, fill=Y)

    l.pack(side=LEFT,fill=BOTH, expand=1)

    l_frame.pack()

    Button(right_drbox, text="Refresh List", command = refList).pack(pady=2)
    Label(right_drbox, text = "").pack(pady=1)

    Label(right_drbox, text = "Output Folder").pack()
    folder = Entry(right_drbox, width=30)
    folder.insert(0,path.get())
    folder.configure(state='readonly')
    folder.pack()

    selFolder = lambda: selectFolder(folder)
    Button(right_drbox, text = 'Change Folder', command = selFolder).pack()

    Label(right_drbox, text = "Output Base Name:").pack()
    o = Entry(right_drbox, width=30)
    o.insert(0,outfile.get())
    o.xview(o.index('end')-20)
    o.pack()

    checkOW = Checkbutton(right_drbox, text="Overwrite Existing File", command=c_OW, takefocus=0)
    checkOW.pack()


    kvfinderSet = lambda: runKVFinder(d, s, p, v, surf, lf, dc, l, o, l_out, p_out,inp_file,lig_file,out_file,l_res,l_out2, right_drbox,right_drbox)
    Label(right_drbox, text = "").pack(pady=1)
    Button(right_drbox, text="Run KVFinder", fg='DarkGreen', command = kvfinderSet).pack(pady=4)

    saveBoxSet = lambda: saveNewBox(d, s, p, v, surf, lf, dc, o, p_out,right_drbox)
    Label(right_drbox, text = "").pack(pady=2)
    Button(right_drbox, text="Save Parameters", fg='Blue', command = saveBoxSet).pack(pady=4)


    deleteBox(min_x, max_x, min_y, max_y, min_z, max_z, ang_s, ang2_s)
    cmd.delete('box')

    refreshList(l)
    refreshList(l_basic)
################ End of the DrawBox Frame ##########################


#################Constructs the Config Tab #########################

    Label(config, text = "").pack(pady="4")

    Label(config, text = "Dictionary Filename:").pack(pady="3")
    d = Entry(config, takefocus=0, width=30)
    d.pack(pady="4")

    restoreConfigSet = lambda: restoreConfig(d)

    Label(config, text = "").pack(pady="8")

    Button(config, text="Restore Default Value", command =restoreConfigSet).pack(pady="5")
################## End of the Config Tab ###########################

################ Constructs the Options Frame ##########################

    left_opt = Frame(opt)
    right_opt = Frame(opt)
    left_opt.pack(side='left', fill='both', padx = 40)
    right_opt.pack(side='left', fill='both', padx = 40)

    Label(left_opt, text = "").pack(pady="4")
    Label(left_opt, text = "Step Size ("+unichr(0x212b)+"):").pack(pady="3")
    s = Entry(left_opt, takefocus=0)
    s.pack(pady="4")

    Label(left_opt, text = "Prob In Size ("+unichr(0x212b)+"):").pack(pady="3")
    p = Entry(left_opt, takefocus=0)
    p.pack(pady="4")


    Label(left_opt, text = "Probe Out Size ("+unichr(0x212b)+"):").pack(pady="3")
    p_out = Entry(left_opt, takefocus=0)
    p_out.pack(pady="2")

    Label(left_opt, text = "").pack(pady='1')


    Label(left_opt, text = "").pack(pady="4")

    c_BMSet = lambda: c_BM(DNBButton,checkWP)
    checkBM = Checkbutton(left_opt, text="Probe Out Adjustment", command=c_BMSet, takefocus=0)
    checkBM.pack()

    c_WPSet = lambda: c_WP(min_x, max_x, min_y, max_y, min_z, max_z,ang_s, ang2_s, DNBButton, checkBM)

    checkWP = Checkbutton(left_opt, text="Whole Protein", command=c_WPSet, takefocus=0)
    checkWP.pack()


    Label(left_opt, text = "").pack(pady="4")


    surf = IntVar()

    selectSASSET = lambda: SelectSAS(surf)

    Label(left_opt, text = "Surface:").pack(pady="3")
    sas = Radiobutton(left_opt, text="Solvent Acessible Surface", value=0, takefocus=0, command=selectSASSET)

    selectMSSET = lambda: SelectMS(surf)

    ms = Radiobutton(left_opt, text="Molecular Surface", value=1, takefocus=0, command=selectMSSET)
    sas.pack()
    ms.pack()

    Label(left_opt, text = "").pack(pady="5")


    checkSR = Checkbutton(left_opt, text="Step Redimensioning", command=c_SR, takefocus=0)
    checkSR.pack()

    Label(right_opt, text = "").pack(pady="4")


    Label(right_opt, text = "Minimum Cavity Volume ("+unichr(0x212b)+unichr(0x00b3)+"):").pack(pady="3")
    v = Entry(right_opt, takefocus=0)
    v.pack(pady="4")

    Label(right_opt, text = "").pack(pady=6)
    c_LASET = lambda: c_LA(lf)
    checkLA = Checkbutton(right_opt, text="Ligand Adjustment", command=c_LASET, takefocus=0)
    checkLA.pack(pady=6)


    Label(right_opt, text = "").pack(pady="3")

    Label(right_opt, text = "Ligand Filename:").pack()

    opt_lf = Frame(right_opt)
    sb_opt_lf = Scrollbar(opt_lf, orient="vertical")
    lf = Listbox(opt_lf, height = 4, width=30, yscrollcommand=sb_opt_lf.set, exportselection=0)

    sb_opt_lf.config(command=lf.yview)
    sb_opt_lf.pack(side=RIGHT, fill=Y)

    lf.pack(side=LEFT,fill=BOTH, expand=1)

    opt_lf.pack()


    selFile = lambda: selectFile(lf)
    refLigFile = lambda: refreshList(lf)

    ligframe = Frame(right_opt)
    Button(ligframe, text = 'Refresh List', command = refLigFile).pack(side="left")
    Button(ligframe, text = 'Another File', command = selFile).pack(side="left")

    ligframe.pack(pady='4')
    Label(right_opt, text = "").pack(pady="4")
    Label(right_opt, text = "Distance Cutoff:").pack(pady="3")
    dc = Entry(right_opt, takefocus=0)
    dc.pack(pady="4")

    Label(right_opt, text = "").pack(pady="4")

    restoreOptionSet = lambda: restoreOptions(s, p, v, sas, ms, checkLA, checkSR, lf, dc, o, p_out, checkOW, checkWP, checkBM, exp, surf)

    Button(right_opt, text="Restore Default Values", command =restoreOptionSet).pack(pady="5")

################ End of the Options Frame ##########################


###################### About Frame #######################################
    Label(about, text = "").pack(pady=10)
    Label(about, text = "The KVFinder software was developed by:").pack(pady=8)
    Label(about, text = "Felipe Augusto Nunes Ferraz").pack()
    Label(about, text = "Saulo Henrique Pires de Oliveira").pack()
    Label(about, text = "Rodrigo Vargas Honorato").pack()
    X_name = "Jos" + unichr(0x00e9) + " Xavier Neto"
    Label(about, text = X_name).pack()

    T_name = "Tiago Jos"+ unichr(0x00e9) +" Paschoal Sobreira"
    Label(about, text = T_name).pack()
    Label(about, text = "Paulo Sergio Lopes de Oliveira").pack()
    Label(about, text = "").pack(pady=10)

    Label(about, text = "National Center of Energy and Material Research - CNPEM\nNational Laboratory of Biosciences - LNBio").pack()
    Label(about, text = "Campinas, Brazil - P.O. Box 6192 - CEP 13083-970, Campinas - SP").pack()


    Label(about, text = "").pack(pady=10)
    Label(about, text = "Contact: paulo.oliveira@lnbio.cnpem.br").pack()
    Label(about, text = "KVFinder website: http://lnbio.cnpem.br/facilities/bioinformatics/software-2/").pack()



#######################End of the about Frame ############################

    restore(d, s, p, v, sas, ms, checkLA, checkSR, lf, dc, o, p_out, checkOW, checkWP, checkBM, exp, s_basic, p_basic, p_out_basic, surf, v_basic)

    for child in mainframe.winfo_children():
        child.grid_configure(padx=5, pady=5, sticky='nwse')
    notebook_adv.setnaturalsize()
    notebook.setnaturalsize()

    root.mainloop();


################ Functions related to the check boxes ###################

def SelectSAS(surf):
    surf.set(0)

def SelectMS(surf):
    surf.set(1)

def c_SR():
    global StepRed
    if StepRed.get() == 0:
        StepRed.set(1)
    else:
        StepRed.set(0)

def c_WP(min_x, max_x, min_y, max_y, min_z, max_z,ang_s, ang2_s, DNBButton, checkWP):
    global WholeProt
    if WholeProt.get() == 0:
        deleteBox(min_x, max_x, min_y, max_y, min_z, max_z,ang_s, ang2_s)
        DNBButton.configure(state='disabled')
        WholeProt.set(1)
        BoxMode.set(0)
        checkWP.deselect()
    else:
        DNBButton.configure(state='normal')
        WholeProt.set(0)

def c_LA(lf):
    global ligandAd
    if ligandAd.get() == 0:
        ligandAd.set(1)
        refreshList(lf)
    else:
        ligandAd.set(0)

def c_OW():
    global OverW
    if OverW.get() == 0:
        OverW.set(1)
    else:
        OverW.set(0)

def c_BM(DNBButton,checkWP):
    global BoxMode
    if BoxMode.get() == 0:
        BoxMode.set(1)
        WholeProt.set(0)
        DNBButton.configure(state='normal')
        checkWP.deselect()
    else:
        BoxMode.set(0)

def c_CD():
    global ChargeDisplay
    control = 0
    if ChargeDisplay.get() == 0:
        ChargeDisplay.set(1)

        for item in cmd.get_names("all"):
            if item == 'cavities.KVFinder.output':
                control = 1
        if control == 0:
            return

        cmd.disable(last_object.get())
        cmd.do("spectrum pc, red_white_blue, cavities.KVFinder.output, minimum=-3, maximum=3")
        cmd.enable(last_object.get())
    else:
        ChargeDisplay.set(0)

        for item in cmd.get_names("all"):
            if item == 'cavities.KVFinder.output':
                control = 1
        if control == 0:
            return

        cmd.disable(last_object.get())
        cmd.color('blue','cavities.KVFinder.output')
        cmd.enable(last_object.get())
