# -*- coding: utf-8 -*-

__author__ = "Konstantinos Vlachas"
__email__ = "vlachask@ibk.baug.ethz.ch"


import numpy as np


class Newmark:
    def __init__(self, Assembly, Parameters=None):
        self.Parameters = Parameters
        self.Assembly = Assembly

    def simulation(self, MODEL=None, Parameters=None):
        if MODEL is None:
            MODEL = self.Assembly

        dt = MODEL.Excitation.dt

        a1, a2, a3, a4, a5, a6 = Newmark.SolverParams(dt, Parameters)

        # Obtain system matrices
        Kstiff=MODEL.K
        M=MODEL.M
        fint = MODEL.fint
        C=MODEL.C

        a1M, a4C = a1 * M, a4 * C

        # Initialize
        Displacements = MODEL.Displacements
        Velocities = MODEL.Velocities
        Accelerations = MODEL.Accelerations
        HystereticCurvesR = MODEL.BoucWenModel.Hist["HystereticCurvesR"]
        HystereticCurvesU = MODEL.BoucWenModel.Hist["HystereticCurvesU"]

        v0 = Velocities[:, 0]
        u0 = Displacements[:, 0]

        # Extract external forcing
        fnorm = MODEL.fnorm
        Fmatrix = MODEL.RHS

        # Compute initial conditions
        f0 = Fmatrix[:, 0]
        a0 = np.linalg.inv(M).dot((f0 - Kstiff.dot(u0) - C.dot(v0)))

        uk = u0
        uik = uk
        vk = v0
        ak = a0

        tol, maxit = 1e-03, 20
        # Start iterations
        for i in range(np.shape(Fmatrix)[1]-1):
            fi = Fmatrix[:, i+1]

            va1 = a2 * vk + a3 * ak
            va2 = a5 * vk + a6 * ak

            ui = uk
            DispsU = ui
            K, Mass, f, BWhistTemp = MODEL.getSystemMatrices(DispsU, Assemble=False)
            Kstiff, M, fint, _, bcdofs = MODEL.applyBCs(K,Mass,f)

            # Compute residual
            R = M.dot(-va1) + C.dot(va2) + fint[:, 0] - fi

            Rnorm = np.linalg.norm(R)
            nit = 0
            uik = ui - uk

            while (Rnorm > tol * fnorm) and (nit <= maxit):
                Keff = a1M + a4C + Kstiff
                # Correct displacement based on residual
                du = np.linalg.inv(Keff).dot(-R)
                ui = ui + du
                DispsU = ui
                uik = ui - uk

                # Update system matrices due to nonlinearity
                K, Mass, f, BWhistTemp = MODEL.getSystemMatrices(DispsU, Assemble=False)
                Kstiff, M, fint, _, bcdofs = MODEL.applyBCs(K,Mass,f)

                # Re-compute residual
                R = M.dot(a1 * uik - va1) + C.dot(a4 * uik + va2) + fint[:, 0] - fi
                Rnorm = np.linalg.norm(R)

                nit = nit + 1

            if nit <= maxit:
                MODEL.BoucWenModel.Hist = BWhistTemp
                if i % 2000 == 0:
                    print(
                        "Step %i converged after %i iterations, residual %f \n"
                        % (i, nit, Rnorm / fnorm)
                    )
            else:
                print(
                    "Step %i did not converge after %i iterations, residual %f \n"
                    % (i, nit, Rnorm / fnorm)
                )
                
                break

            Displacements[:, i + 1] = DispsU
            uk = ui
            ak = a1 * uik - va1
            vk = a4 * uik + va2
            
            vk[bcdofs] = 0
            ak[bcdofs] = 0
            
            Velocities[:, i + 1] = vk
            Accelerations[:, i + 1] = ak
            HystereticCurvesU[:, i + 1] = MODEL.BoucWenModel.Hist["HistU"]  # HistDict["HistZeta"]
            HystereticCurvesR[:, i + 1] = MODEL.BoucWenModel.Hist["HistR"]

        Results = dict()
        Results["Displacements"] = Displacements
        Results["Velocities"] = Velocities
        Results["Accelerations"] = Accelerations
        Results["HystereticR"] = HystereticCurvesR
        Results["HystereticU"] = HystereticCurvesU
        originaldofs=6*(MODEL.Nodes.shape[0] - MODEL.nl_link_elements.shape[0])        
        Results["OriginalDisps"] = Displacements[:originaldofs,:]

        return Results

    def SolverParams(dt, Parameters=None):
        if Parameters is None:
            gamma, beta=1/2, 1/4
        else:
            gamma=Parameters.gamma
            beta=Parameters.beta

        # Newmark parameters
        a1 = 1 / (beta * (dt ** 2))
        a2 = 1 / (beta * dt)
        a3 = (1 - 2 * beta) / (2 * beta)
        a4 = gamma / (beta * dt)
        a5 = 1 - gamma / beta
        a6 = (1 - gamma / (2 * beta)) * dt

        return a1, a2, a3, a4, a5, a6