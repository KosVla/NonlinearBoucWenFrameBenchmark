# -*- coding: utf-8 -*-

__author__ = "Konstantinos Vlachas"
__email__ = "vlachask@ibk.baug.ethz.ch"


import numpy as np
from matplotlib import pyplot as plt


 
def Plot_TH_Response(Field, t, ndof_plot=None, title="Time History Response"):
    if ndof_plot is None:
        ndof_plot = np.argmax(np.amax(abs(Field), axis=1))
    plt.plot(t, Field[ndof_plot, :], color="blue", label=title)

def Plot_Hysteretic_Curve(HistU, HistR, ndof_plot=None):    
    if ndof_plot is None:
        ndof_plot = np.argmax(np.amax(abs(HistR), axis=1))  
    
    plt.plot(HistU[ndof_plot, :], HistR[ndof_plot, :], color="blue", label="Hysteretic Curve")  


def Plot_Phase_Space(U,V, ndof_plot=None):
    if ndof_plot is None:
        ndof_plot = np.argmax(np.amax(abs(V), axis=1))  
        
    plt.plot(U[ndof_plot, :], V[ndof_plot, :], color="blue", label="Phase Space")