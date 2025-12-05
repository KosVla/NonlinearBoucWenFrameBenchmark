# -*- coding: utf-8 -*-

__author__ = "Konstantinos Vlachas"
__email__ = "vlachask@ibk.baug.ethz.ch"


import numpy as np

# import Assembly

class BoucWenModel:
    def __init__(self,bw_k=[2.0e8,2.0e8], bwa=0.25, Beta=3.,
        Gamma=2., Alpha=1.0, N=1., deltav=0., deltan=0.,
        integration_method='RK4'):
        
        self.Models= BoucWenModel.set_BW_Models(bw_k, bwa, Beta,
        Gamma, Alpha, N, deltav, deltan)
        self.integration_method=integration_method

    def set_BW_Models(bw_k, bwa, Beta,Gamma, Alpha, N, deltav,
                      deltan):

        if isinstance(bwa, list):
            bwas=bwa
        else: 
            bwas = [bwa for x in bw_k]
            
        Alphas = [Alpha for x in bw_k]
        Betas = [Beta for x in bw_k]
        Gammas = [Gamma for x in bw_k]
        Ns = [N for x in bw_k]
        deltavs = [deltav for x in bw_k]
        deltans = [deltan for x in bw_k]

        if isinstance(bw_k, list):
            Models = np.array(list(zip(bwas, bw_k, Alphas, Betas, Gammas, Ns,
                                        deltavs, deltans)))
        else:
            Models = [bwas, bw_k, Alpha, Beta, Gamma, N, deltav, deltan]

        return Models
    
    def change_BWks(self, bw_k, nmaterial=None):

        if nmaterial is None:
            self.Models[:,1]=bw_k
        else:
            self.Models[nmaterial,1]=bw_k

        self.Models

    
    def AssembleProperties(self, nl_links):
        Properties = np.zeros(shape=(nl_links.shape[0] * 6, 10))
        count = 0

        for e in range(nl_links.shape[0]):
            link_nodes = nl_links[e, :2]
            link_dofs = BoucWenModel.elementglobaldofs(link_nodes)
            link_type = nl_links[e, -1]
            
            for d in range(6):
                Properties[count, 0:2] = [link_dofs[d], link_dofs[d + 6]]
                nl_links_properties = self.Models[link_type,:]
                Properties[count, 2:] = nl_links_properties

                count += 1

        self.Properties = Properties
    
    def elementglobaldofs(nodes):
        dofs = np.zeros(shape=(len(nodes) * 6))
        dofs[0: len(nodes) * 6: 6] = (nodes - 1) * 6
        dofs[1: len(nodes) * 6: 6] = (nodes - 1) * 6 + 1
        dofs[2: len(nodes) * 6: 6] = (nodes - 1) * 6 + 2
        dofs[3: len(nodes) * 6: 6] = (nodes - 1) * 6 + 3
        dofs[4: len(nodes) * 6: 6] = (nodes - 1) * 6 + 4
        dofs[5: len(nodes) * 6: 6] = (nodes - 1) * 6 + 5
        return dofs
    
    def evaluateBoucWen(self, U, BWhistory):

        Properties=self.Properties
        Ulinks = U[Properties[:, 0].astype(int)] - U[Properties[:, 1].astype(int)]
    
        if self.integration_method=='Euler':
            return self.evaluateBoucWenEuler(Ulinks, BWhistory)
        elif self.integration_method=='RK4':
            return self.evaluateBoucWenRK4(Ulinks, BWhistory)
        elif self.integration_method=='RK2':
            return self.evaluateBoucWenRK(Ulinks, BWhistory)
        else:
            raise ValueError(f"Unknown integration method: {self.integration_method}")
        
    def evaluateBoucWenEuler(self, U, historyBW):

        histR = historyBW["HistR"]
        histU = historyBW["HistU"]
        histZeta = historyBW["HistZeta"]
        histE = historyBW["HistE"]
        
        deltaU = U - histU

        Properties=self.Properties
        bwa, bw_k, Alpha, Beta, Gamma, N, deltav, deltan=Properties[:,2],Properties[:,3],Properties[:,4],Properties[:,5],Properties[:,6],Properties[:,7],Properties[:,8],Properties[:,9]

        # Strength deterioration
        Energy = histE
        ve = 1 + deltav * Energy
        # Stiffness degradation
        ne = 1 + deltan * Energy

        zetaprior = histZeta
        kbw = (
            Alpha
            - ve
            * (Beta * np.sign(zetaprior) * np.sign(deltaU) + Gamma)
            * np.absolute(zetaprior) ** (N)
        ) / ne

        # Force terms
        LinearTerm = bwa * bw_k * deltaU
        dzeta = kbw * deltaU
        HystereticTerm = (1 - bwa) * bw_k * dzeta

        # Restoring force
        r = histR + HystereticTerm + LinearTerm
        k = bwa * bw_k + (1 - bwa) * bw_k * kbw

        # Update History
        HistR = r
        HistU = histU + deltaU
        HistZeta = zetaprior + dzeta
        #Correction in energy term
        HistE = histE + (zetaprior + HistZeta) / 2.0 * deltaU 

        bwhistory = dict()
        bwhistory["HistR"] = HistR
        bwhistory["HistU"] = HistU
        bwhistory["HistZeta"] = HistZeta
        bwhistory["HistE"] = HistE

        return r, k, bwhistory
    
    def evaluateBoucWenRK(self, Unode, historyBW):

        histR = historyBW["HistR"]
        histU = historyBW["HistU"]
        histZeta = historyBW["HistZeta"]
        histE = historyBW["HistE"]

        Properties=self.Properties
        bwa, bw_k, Alpha, Beta, Gamma, N, deltav, deltan=Properties[:,2],Properties[:,3],Properties[:,4],Properties[:,5],Properties[:,6],Properties[:,7],Properties[:,8],Properties[:,9]

        #Equivalent to Dt
        deltaU = Unode - histU

        # Strength deterioration
        Energy = histE
        ve = 1 + deltav * Energy
        # Stiffness degradation
        ne = 1 + deltan * Energy

        zetaprior = histZeta
        kbw = self.compute_kbw(zetaprior, deltaU, Alpha, Beta, Gamma, N, ve, ne)
        RK_k1 = kbw * deltaU

        zetaprior_k2=zetaprior+0.5*RK_k1
        kbw_k2 = self.compute_kbw(zetaprior_k2, deltaU, Alpha, Beta, Gamma, N, ve, ne)

        RK_k2=kbw_k2*deltaU


        # Force terms
        LinearTerm = bwa * bw_k * deltaU
        dzeta = RK_k2
        HystereticTerm = (1 - bwa) * bw_k * dzeta

        # Restoring force
        r = histR + HystereticTerm + LinearTerm
        # Evaluate kbw at final state
        zeta_final = zetaprior + dzeta
        kbw_final = self.compute_kbw(zeta_final, deltaU,Alpha, Beta, Gamma, N, ve, ne)
        k = bwa * bw_k + (1 - bwa) * bw_k * kbw_final

        # Update History
        HistR = r
        HistU = histU + deltaU
        HistZeta = zeta_final
        #Correction in energy term
        HistE = histE + (zetaprior + HistZeta) / 2.0 * deltaU 

        bwhistory = dict()
        bwhistory["HistR"] = HistR
        bwhistory["HistU"] = HistU
        bwhistory["HistZeta"] = HistZeta
        bwhistory["HistE"] = HistE

        return r, k, bwhistory
    
    def evaluateBoucWenRK4(self, Unode, historyBW):
        """
        RK4 integration of Bouc-Wen model
        4th-order accurate: O(ΔU^4)
        """
        Properties=self.Properties
        bwa, bw_k, Alpha, Beta, Gamma, N, deltav, deltan=Properties[:,2],Properties[:,3],Properties[:,4],Properties[:,5],Properties[:,6],Properties[:,7],Properties[:,8],Properties[:,9]

        histR = historyBW["HistR"]
        histU = historyBW["HistU"]
        histZeta = historyBW["HistZeta"]
        histE = historyBW["HistE"]

        #Equivalent to Dt
        deltaU = Unode - histU

        # Strength deterioration
        Energy = histE
        ve = 1 + deltav * Energy
        # Stiffness degradation
        ne = 1 + deltan * Energy

        zetaprior = histZeta
        # RK4 stages
        k1 = self.compute_kbw(zetaprior, deltaU,  Alpha, Beta, Gamma, N, ve, ne) * deltaU

        k2 = self.compute_kbw(zetaprior + 0.5*k1, deltaU, Alpha, Beta, Gamma, N, ve, ne) * deltaU

        k3 = self.compute_kbw(zetaprior + 0.5*k2, deltaU, Alpha, Beta, Gamma, N, ve, ne) * deltaU

        k4 = self.compute_kbw(zetaprior + k3, deltaU, Alpha, Beta, Gamma, N, ve, ne) * deltaU

        # Weighted average (RK4 formula)
        dzeta = (k1 + 2*k2 + 2*k3 + k4) / 6.0

        # Force terms
        LinearTerm = bwa * bw_k * deltaU
        HystereticTerm = (1 - bwa) * bw_k * dzeta

        # Tangent stiffness - evaluate at final state for consistency
        r = histR + HystereticTerm + LinearTerm
        zeta_new = zetaprior + dzeta
        kbw_final = self.compute_kbw(zeta_new, deltaU, Alpha, Beta, Gamma, N, ve, ne)
        k = bwa * bw_k + (1 - bwa) * bw_k * kbw_final

        # Update History
        HistR = r
        HistU = histU + deltaU
        HistZeta = zeta_new
        #Correction in energy term
        HistE = histE + (zetaprior + HistZeta) / 2.0 * deltaU 

        bwhistory = dict()
        bwhistory["HistR"] = HistR
        bwhistory["HistU"] = HistU
        bwhistory["HistZeta"] = HistZeta
        bwhistory["HistE"] = HistE

        return r, k, bwhistory
    
    def compute_kbw(self, zeta, delta_u, Alpha, Beta, Gamma, N, ve, ne):
        
        return (Alpha- ve * (Beta * np.sign(zeta) * np.sign(delta_u)+ Gamma)* np.absolute(zeta) ** N) / ne



        
    
