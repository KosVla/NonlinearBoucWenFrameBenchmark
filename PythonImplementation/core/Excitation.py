# -*- coding: utf-8 -*-

__author__ = "Konstantinos Vlachas"
__email__ = "vlachask@ibk.baug.ethz.ch"

import numpy as np

class Excitation:
    def __init__(
        self,
        fintegration=1000,
        angle=np.pi / 4,
        Amp=10,
        type="Kobe",
    ):

        # integration sampling frequency.
        self.dt = 1 / fintegration
        # integration time step.
        self.angle = angle
        self.Amp = Amp

        if type == "Kobe":
            Excite = np.loadtxt("InputFiles/KobeAccelNoScaling.txt", delimiter=",")
            Excitation = Amp*Excite[1,:]
            tsample = Excite[0,:]
        elif type == "ElCentro":
            Excite = np.loadtxt("InputFiles/ElCentroAccelNoScaling.txt", delimiter=",")
            Excitation = Excite[1,:]
            tsample = Amp*Excite[0,:]
        elif type == "Morgan":
            Excite = np.loadtxt("InputFiles/MorganAccelNoScaling.txt", delimiter=",")
            Excitation = Amp*Excite[1,:]
            tsample = Excite[0,:]

        xq = np.arange(0, tsample[-1], 1/fintegration)        

        if tsample[1]<1/fintegration:
            print('Warning: Integration timestep is larger than timestep of recordings')

        SynthesizedAccelerogram = np.interp(xq, tsample, Excitation)
        self.SynthesizedAccelerogram = np.reshape(SynthesizedAccelerogram,(np.shape(SynthesizedAccelerogram)[0],1))
        self.time = xq