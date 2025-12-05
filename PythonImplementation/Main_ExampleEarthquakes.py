# -*- coding: utf-8 -*-
__author__ = "Konstantinos Vlachas"
__email__ = "vlachask@ibk.baug.ethz.ch"

################################
# Main file
################################

################################
# Import necessary libraries
################################
import numpy as np
from matplotlib import pyplot as plt
from InputFiles.InputFile import *
from core.Utils import *
from core.Assembly import *
from core.Excitation import *
from core.Newmark import *
from core.BoucWenModel import *

################################
# Define Model Parameters
################################
InputFile='BeamsExample'
Mat = Material([210e9, 210e9, 210e9, 210e9])
Forcing = Excitation(fintegration=1000, Amp=1e6, angle=np.pi/4, type="Morgan")
bw_ks = [1.0, 1.0, 1.0, 1.0] #Basement, Columns, First, Second
BW = BoucWenModel([x*2.e8 for x in bw_ks])

################################
# Instatiate model
################################
MODEL = ModelAssembly(Mat, Forcing, BW, np.ones((1,1)),InputFile)

################################
# Integrate Model
################################
Solver = Newmark(MODEL)

Results = Solver.simulation()

################################
# Plot Results
################################

#Plot_TH_Response(Results["Disps"], SimulationParameters.time, title="Displacement Response")
Plot_TH_Response(Results["Accelerations"], Forcing.time, title="Acceleration Response")
Plot_TH_Response(Results["OriginalDisps"], Forcing.time, title="Displacement Response")
Plot_Hysteretic_Curve(Results["HystereticU"], Results["HystereticR"])
Plot_Phase_Space(Results["Displacements"], Results["Velocities"])

print(ARK)