# -*- coding: utf-8 -*-

__author__ = "Konstantinos Vlachas"
__email__ = "vlachask@ibk.baug.ethz.ch"


import numpy as np
import scipy.io
import random
from matplotlib import pyplot as plt
import scipy as sp
import CodeDamaged.POD_SVDtruncBasis as BAT
from scipy.signal import butter,filtfilt
from scipy import interpolate


class Utilities:
    def __init__(self, Parameters=None):
        self.parameters = Parameters

    @staticmethod
    def GetFOMresponse(ROMresponse, Basis):

        Ref = np.zeros(shape=(Basis.shape[0], ROMresponse.shape[1]))

        for i in range(ROMresponse.shape[1]):
            Ref[:, i] = Basis.dot(ROMresponse[:, i])

        return Ref
 

    @staticmethod
    def PoluteSignal(CleanSignals, delta=None):
        #CleanSignals is (s x n):
        #   n are the timesteps, s are the dofs
        CleanSignals = CleanSignals.T

        mean = np.mean(CleanSignals, axis=0)
        std = np.sqrt(np.sum((CleanSignals-mean)**2/CleanSignals.shape[0], axis=0))

        noisy_signal=CleanSignals+delta*std*np.random.normal(size=CleanSignals.shape)
        noisy_signal=noisy_signal.squeeze()

        return noisy_signal.T

    @staticmethod
    def PlotResponse(ResultsH, ResultsD=None):
        Uref = ResultsH["Disps"]        
        ndof_plot = np.argmax(np.amax(abs(Uref), axis=1))

        # Plot tip
        plt.plot(Uref[ndof_plot, :], color="blue", label='Reference (FOM) response')
        if ResultsD is not None:
            U = ResultsD["Disps"]
            plt.plot(U[ndof_plot, :], color="orange", linestyle="dashed", label="ROM approximation")
        
        plt.xlabel("Time increments")
        plt.ylabel("Displacement Response")
        plt.legend(loc='lower right')
        plt.tight_layout()
        plt.show()

    @staticmethod
    def PlotResponseA(ResultsH, ResultsD=None, Basis=None):
        Uref = ResultsH["AccelerationsA"]        
        ndof_plot = np.argmax(np.amax(abs(Uref), axis=1))

        # Plot tip
        plt.plot(Uref[ndof_plot, :], color="blue", label='Reference (FOM) response')
        if ResultsD is not None:
            U = ResultsD["AccelerationsA"]
            U = Utilities.GetFOMresponse(U, Basis)
            plt.plot(U[ndof_plot, :], color="orange", linestyle="dashed", label="ROM approximation")
        
        plt.xlabel("Time increments")
        plt.ylabel("Acceleration Response")
        plt.legend(loc='lower right')
        plt.tight_layout()
        plt.show()

    @staticmethod
    def PlotHysteretic(ResultsH, ResultsD=None, ndof=None):
        # Plot Hysteresis curve
        Rh = ResultsH["HystereticParamsR"]
        Uh = ResultsH["HystereticParamsU"]

        if ndof is None:
            # Pick dof
            ndof_plot = np.argmax(np.amax(abs(Uh), axis=1))
        else:
            ndof_plot = ndof

        # Plot
        plt.plot(Uh[ndof_plot, :], Rh[ndof_plot, :], color="blue")
        if ResultsD is not None:
            R = ResultsD["HystereticParamsR"]      
            U = ResultsD["HystereticParamsU"]
            plt.plot(U[ndof_plot, :], R[ndof_plot, :], color="red", linestyle="dashed")
    
        plt.xlabel("Incremental Displacements")
        plt.ylabel("Hysteretic Forcing")
        plt.tight_layout()
        plt.show()

    @staticmethod
    def ClosestBasis(Sample,Parameters,BasesPool,k=1,distances='euclid'):
        #Compute distances
        if distances == 'euclid':
            Sample = Sample.reshape((1,len(Sample)))
            dist = (Parameters - Sample)**2
            dist = np.sum(dist, axis=1)
            dist = np.sqrt(dist)

        elif distances == 'mahalanobis':
            # Covariance matrix
            covariance = np.cov(Parameters, rowvar=True)

            # Covariance matrix power of -1
            # covariance_pm1 = np.linalg.matrix_power(covariance, -1)
            inv_covmat = sp.linalg.inv(covariance)
            # Center point
            centerpoint = np.mean(Parameters, axis=1)     
            x_minus_mu = (
                Sample- np.tile(centerpoint.T, (np.shape(Sample)[1], 1)).T
            )
            left_term = np.dot(x_minus_mu.T, inv_covmat)
            mahal = np.dot(left_term, x_minus_mu)
            md = np.sqrt(mahal.diagonal())
            

        #Find closest index
        index = np.argsort(dist)

        #Select and return basis
        if k==1:
            LocalPODbasis = BasesPool[index[0],:,:]
        else:
            Vlocals=BasesPool[index[0],:,:]
            for i in range(k-1):
                Basis = BasesPool[index[i+1],:,:]
                Vlocals = np.hstack((Vlocals,Basis))
            LocalPODbasis = BAT.PODbasis(Vlocals, np.shape(BasesPool)[2])
        
        return LocalPODbasis
        


    @staticmethod
    def ErrorROMFOM(ResultsRef, ResultsVal, Basis):

        UrefR = ResultsRef["Disps"]
        Uval = ResultsVal["Disps"]
        VrefR = ResultsRef["Velocities"]
        Vval = ResultsVal["Velocities"]
        ArefR = ResultsRef["AccelerationsA"]
        Aval = ResultsVal["AccelerationsA"]
        
        Rref = ResultsRef["HystereticParamsR"]
        Rval = ResultsVal["HystereticParamsR"]
        ErrorR = np.linalg.norm(Rref - Rval) / np.linalg.norm(Rval)
        ndof_max = np.argmax(np.amax(abs(Rval), axis=1))
        ErrorRmax = np.linalg.norm(
            Rref[ndof_max, :] - Rval[ndof_max,:]
        ) / np.linalg.norm(Rval[ndof_max, :])

        Vref = np.zeros(shape=(Vval.shape[0], Vval.shape[1]))
        Aref = np.zeros(shape=(Vval.shape[0], Vval.shape[1]))

        if UrefR.shape[0] != Uval.shape[0]:
            Uref = np.zeros(shape=(Uval.shape[0], Uval.shape[1]))
        else:
            Uref = UrefR

        for i in range(Uval.shape[1]):
            Vref[:, i] = Basis.dot(VrefR[:, i])
            Aref[:, i] = Basis.dot(ArefR[:, i])
            if UrefR.shape[0] != Uval.shape[0]:
                Uref[:, i] = Basis.dot(UrefR[:, i])

        ErrorU = np.linalg.norm(Uref - Uval) / np.linalg.norm(Uval)
        ErrorV = np.linalg.norm(Vref - Vval) / np.linalg.norm(Vval)
        ErrorA = np.linalg.norm(Aref - Aval) / np.linalg.norm(Aval)

        ndof_max = np.argmax(np.amax(abs(Uval), axis=1))
        ErrorUmax = np.linalg.norm(
            Uref[ndof_max, :] - Uval[ndof_max, :]
        ) / np.linalg.norm(Uval[ndof_max, :])
        # ndof_max = np.argmax(np.amax(abs(Vval), axis=1))
        # ErrorVmax = np.linalg.norm(
        #    Vref[ndof_max, :] - Vval[ndof_max, :]
        # ) / np.linalg.norm(Vval[ndof_max, :])
        ndof_max = np.argmax(np.amax(abs(Aval), axis=1))
        ErrorAmax = np.linalg.norm(
            Aref[ndof_max, :] - Aval[ndof_max, :]
        ) / np.linalg.norm(Aval[ndof_max, :])

        Errors = {}
        Errors["ErrorU"] = ErrorU
        Errors["ErrorV"] = ErrorV
        Errors["ErrorA"] = ErrorA
        Errors["ErrorR"] = ErrorR
        Errors["ErrorUmax"] = ErrorUmax
        # Errors["ErrorVmax"] = ErrorVmax
        Errors["ErrorAmax"] = ErrorAmax
        Errors["ErrorRmax"] = ErrorRmax
        Errors["SpeedUp"] = ResultsVal["Time"] / ResultsRef["Time"]
        #Errors["VelocitiesFull"]=Vref

        return Errors
    
    @staticmethod
    def ErrorROMFOMA(ResultsRef, ResultsVal, Basis):

        ArefR = ResultsRef["AccelerationsA"]
        Aval = ResultsVal["AccelerationsA"]
        
        Aref = np.zeros(shape=(Aval.shape[0], Aval.shape[1]))
        for i in range(Aval.shape[1]):
            Aref[:, i] = Basis.dot(ArefR[:, i])

        ndofs= np.argsort(np.max(abs(Aval), axis=1), axis=0)
        compute = ndofs[-10:]

        ErrorAmax = np.linalg.norm(
            Aref[compute, :] - Aval[compute, :]
        ) / np.linalg.norm(Aval[compute, :])

        return ErrorAmax

    @staticmethod
    def ErrorFOM(ResultsL, ResultsT, option=0):
        Uval = ResultsT["Disps"]
        Ur = ResultsL["Disps"]
        Vval = ResultsT["Velocities"]
        Vr = ResultsL["Velocities"]
        Aval = ResultsT["AccelerationsA"]
        Ar = ResultsL["AccelerationsA"]

        if option > 0:
            Rval = ResultsT["HystereticParamsR"]
            Rr = ResultsL["HystereticParamsR"]
            ErrorRh = np.linalg.norm(Rval - Rr) / np.linalg.norm(Rval)
        else:
            ErrorRh = 0

        ErrorU = np.linalg.norm(Uval - Ur) / np.linalg.norm(Uval)
        ErrorV = np.linalg.norm(Vval - Vr) / np.linalg.norm(Vval)
        ErrorA = np.linalg.norm(Aval - Ar) / np.linalg.norm(Aval)

        Errors = {}
        Errors["ErrorU"] = ErrorU
        Errors["ErrorV"] = ErrorV
        Errors["ErrorA"] = ErrorA
        Errors["ErrorRh"] = ErrorRh

        return Errors

    @staticmethod
    def Multisine(fmin, fmax, fsample, fstep=1, Ncomponents=1000):
                   
        vals = np.arange(fmin, fmax+1, fstep)

        t = np.arange(0, Ncomponents*1/fsample, 1 / fsample)
        x = np.zeros((len(t),))

        for k in vals:
            phiinit = random.uniform(0, 2 * np.pi)
            x += np.sin(2 * np.pi * k * t + phiinit)

        #x = x / np.amax(abs(x)) * random.uniform(0.2, 1)
        x = x / np.amax(abs(x))

        return x

    @staticmethod
    def fftnoise(f, seedex=30):
        f = np.array(f, dtype="complex")
        Np = (len(f) - 1) // 2
        np.random.seed(seedex)
        #random.seed(seedex)
        #phases = random.random()*2*np.pi
        phases = np.random.rand(Np) * 2 * np.pi
        phases = np.cos(phases) + 1j * np.sin(phases)
        f[1 : Np + 1] *= phases
        f[-1 : -1 - Np : -1] = np.conj(f[1 : Np + 1])

        return np.fft.ifft(f).real

    @staticmethod
    def WhiteNoise(min_freq, max_freq, fsample, samples=1024, randomseed=30):
        fsamples=(max_freq-min_freq)*10
        freqs = np.abs(np.fft.fftfreq(fsamples, 1 / fsample))
        f = np.zeros(samples)
        idx = np.where(np.logical_and(freqs >= min_freq, freqs <= max_freq))[0]
        f[idx] = 1
        y = Utilities.fftnoise(f, randomseed)
        

        return y

    @staticmethod
    def butter_lowpass_filter(data, cutoff, fs, order):
        nyq = fs/2
        normal_cutoff = cutoff / nyq
        # Get the filter coefficients 
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        y = filtfilt(b, a, data)
        
        return y


class SimulationParameters:
    def __init__(
        self,
        fs=100,
        angle=np.pi / 4,
        AmpF=10,
        type="ExampleB",
        fminimum=3,
        fmax=10,
        fstep=1,
        Npoints=1000,
        seedrnd=30,
    ):

        # integration sampling frequency.
        self.dt = 1 / fs
        # integration time step.
        self.angle = angle

        self.lowpass = False

        if type == "ExampleA":
            #Excite = scipy.io.loadmat("ExampleA.mat")
            #Rth = np.array(Excite.get("ExciteA"))
            Excite = np.loadtxt("ExampleA.txt", delimiter=",")
            Rth = np.reshape(Excite[1250:],(1798,1))
            #Rth = np.reshape(Excite,(3048,1))
        elif type == "ExampleB":
            Excite = scipy.io.loadmat("ExampleB.mat")
            R = np.array(Excite.get("ExciteB"))
            Rth = R[0:2000]
            #f = interpolate.interp1d(np.reshape(np.arange(0,15,1e-2),(1500,)), R.T)
            #xnew = np.arange(0, 14.5, 5e-3)
            #Rth = f(xnew) 
            #Rth=Rth.T
            
        elif type == "multisine":
            # Multisine constructor
            Rth = Utilities.Multisine(fminimum, fmax, fs, Npoints, fstep)
        elif type == "whitenoise":
            # White Noise constructor
            Rth = Utilities.WhiteNoise(fminimum, fmax, fs, samples=Npoints)
            Rth = Rth / np.amax(abs(Rth))
        elif type=="manual":
            Rth = Npoints

        self.SynthesizedAccelerogram = AmpF * Rth