# -*- coding: utf-8 -*-

__author__ = "Konstantinos Vlachas"
__email__ = "vlachask@ibk.baug.ethz.ch"


import numpy as np
import numpy.linalg as ln
from InputFiles.InputFile import *
from core.BoucWenModel import *
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigs


class ModelAssembly:

    def __init__(self, Material, Excitation, BoucWenModel, perturb=np.ones((1,1)), InputFile='FloorsExample'):
        AssembleFrame(self, InputFile)
        self.Material = Material
        self.BoucWenModel = BoucWenModel
        self.Excitation = Excitation
        if np.sum(perturb)>1:
            self.perturb=perturb
        else:
            self.perturb=np.ones((self.BeamElements.shape[0],1))
        self.initialize()


    def initialize(self):
        
        self.BoucWenModel.AssembleProperties(self.nl_link_elements)
        nlinks = self.nl_link_elements.shape[0] * 6

        d = dict()
        d["HistR"] = np.zeros(shape=(nlinks))
        d["HistU"] = np.zeros(shape=(nlinks))
        d["HistZeta"] = np.zeros(shape=(nlinks))
        d["HistE"] = np.zeros(shape=(nlinks))

        # Initial conditions
        ndofs = self.n_dofs
        nt = np.shape(self.Excitation.SynthesizedAccelerogram)[0]
        
        DisplacementsU = np.zeros(shape=(ndofs, nt ))
        VelocitiesV = np.zeros(shape=(ndofs, nt ))
        AccelerationsA = np.zeros(shape=(ndofs, nt ))
        d["HystereticCurvesR"] = np.zeros(shape=(nlinks, nt ))
        d["HystereticCurvesU"] = np.zeros(shape=(nlinks, nt ))
        
        self.BoucWenModel.Hist = d

        # Assemble system matrices (stiffness,mass, internal forces)
        Kstiff, Mmass, fint, BWhist = self.getSystemMatrices()
        Kmat, Mmat, Fint, freedofs, bcdofs = self.applyBCs(
            Kstiff,
            Mmass,
            fint,
        )
        #self.BoucWenModel.Hist=BWhist

        # Damping treatment
        if self.ICs["zeta"] is not None:
            Kd=Kstiff[freedofs,:][:,freedofs]
            Md=Mmass[freedofs,:][:,freedofs]

            alpha, beta = self.RayleighDamping(
                Kd, Md,
                self.ICs["zeta"],
                self.ICs["OmegaIndexes"],
            )

        else:
            alpha, beta = 0, 0.05

        Crayleigh = alpha * Mmat + beta * Kmat

        self.K=Kmat
        self.M=Mmat
        self.C=Crayleigh
        self.fint=Fint

        if self.ICs is not None:
            VelocitiesV[:, 0] = self.ICs["Velocities"][:, 0]
            DisplacementsU[:, 0] = self.ICs["Disps"][:, 0]

        self.Displacements = DisplacementsU
        self.Velocities = VelocitiesV
        self.Accelerations = AccelerationsA

        # Assembly of excitation vector
        Zacceler, S = self.AssembleSandZ(Mmass)

        Fmatrix = np.zeros(shape=(ndofs, nt))
        Fmatrix[:, 0: np.shape(Zacceler)[0]] = -S.dot(Zacceler.T)

        self.RHS = Fmatrix
        self.pos = np.argmax(Zacceler)
        self.fnorm = np.linalg.norm(Fmatrix[:, self.pos])

    def getSystemMatrices(self, DispsU=None, HistDict=None, Assemble=True):
        ndof=self.n_dofs
        
        if DispsU is None:
            DispsU = np.zeros(shape=(ndof))

        if Assemble:
            Kstiff, Mmass = np.zeros(
            shape=(ndof, ndof)), np.zeros(shape=(ndof, ndof))
            fint = np.zeros(shape=(ndof, 1))
            # Assembly of beam elements
            self.assembleBeams(DispsU, Kstiff, Mmass, fint)
        else:
            Kstiff, Mmass, _ = self.Kbeams, self.Mbeams, self.fbeams
            fint = Kstiff.dot(DispsU).reshape(-1, 1)
              
        # Assembly of nonlinear links
        if HistDict is None:
            HistDict = self.BoucWenModel.Hist

        Rs, Ks, BWhist = self.BoucWenModel.evaluateBoucWen(DispsU, HistDict)

        DofsLoop = self.BoucWenModel.Properties

        rows = np.concatenate((DofsLoop[:, 0], DofsLoop[:, 1]), axis=None)
        cols = np.zeros(shape=(2 * DofsLoop[:, 0].shape[0]))
        valuesf = np.concatenate((Rs, -Rs), axis=None)
        Fnl = coo_matrix((valuesf, (rows, cols)), shape=(ndof, 1))

        fint = fint + Fnl.toarray()

        rows = np.concatenate(
            (DofsLoop[:, 0], DofsLoop[:, 0],
             DofsLoop[:, 1], DofsLoop[:, 1]), axis=None
        )
        cols = np.concatenate(
            (DofsLoop[:, 0], DofsLoop[:, 1],
             DofsLoop[:, 0], DofsLoop[:, 1]), axis=None
        )
        values = np.concatenate((Ks, -Ks, -Ks, Ks), axis=None)
        Kbw = coo_matrix((values, (rows, cols)), shape=(ndof, ndof))

        Kstiff = Kstiff + Kbw.toarray()
       
        return Kstiff, Mmass, fint, BWhist
    
    def assembleBeams(self, DispsU, Kstiff, Mmass, fint):
        
        beam_cross = self.Material.cross_sections

        for e in range(self.BeamElements.shape[0]):
            el = self.BeamElements[e]
            beam_nodes = el[1:]
            beam_nodes = beam_nodes.astype(int)
            beam_Coords = self.Nodes[beam_nodes - 1, :]

            beam_dofs = self.elementglobaldofs(beam_nodes)
            beam_dofs = beam_dofs.astype(int)
                
            el_type = self.BeamElements[e,0]
            beam_mat_props = self.Material.mat_properties[el_type,:]

            beam_mat_props[0]=beam_mat_props[0]*self.perturb[e]

            Ke, Me, fe = self.assembleBeamElement(
                beam_mat_props, beam_cross, beam_Coords
            )

            Ue = DispsU[beam_dofs]

            Kstiff[np.ix_(beam_dofs, beam_dofs)] = (
                Kstiff[np.ix_(beam_dofs, beam_dofs)] + Ke
            )

            Mmass[np.ix_(beam_dofs, beam_dofs)] = (
                Mmass[np.ix_(beam_dofs, beam_dofs)] + Me
            )

            fint[beam_dofs, 0] = fint[beam_dofs, 0] + Ke.dot(Ue)

        self.Kbeams=Kstiff
        self.Mbeams=Mmass
        self.fbeams=fint

       
    def applyBCs(self, Kstiff, Mmass, fint):
        ndofs = self.n_dofs

        nodal_disps = self.NodalDisplacements
        bc_dofs = []
        bc_dofs_nz = []
        bc_values = []

        for i in range(nodal_disps.shape[0]):
            node = nodal_disps[i, 0]
            dofs = np.arange(6 * (node - 1), 6 * node)
            dofs_act = dofs[np.where(nodal_disps[i, 1:7] != 0)[0]]
            bc_dofs.extend(dofs_act)
            nz = np.where(nodal_disps[i, 7:] != 0)[0]
            if nz.size > 0:
                bc_dofs_nz.extend(dofs[nz])
                bc_values.extend(nodal_disps[i, nz + 7])

        bc_dofs_nzA = np.array(bc_dofs_nz)
        bc_valuesA = np.array(bc_values)

        #if bc_dofs_nzA.size > 0:
        #    fbc = Kstiff[:, bc_dofs_nzA].dot(bc_valuesA)
        #    f = f - fbc

        bcdofs = np.array(bc_dofs, dtype=int)

        Kstiff[bcdofs, :] = 0
        Kstiff[:, bcdofs] = 0
        Mmass[bcdofs, :] = 0
        Mmass[:, bcdofs] = 0
        # f[bcdofs, 0] = 0
        fint[bcdofs, 0] = 0

        if bc_dofs_nzA.size > 0:
            #f[bc_dofs_nzA] = bc_valuesA
            fint[bc_dofs_nzA] = 0

        Kstiff[bcdofs, bcdofs] = 1
        Mmass[bcdofs, bcdofs] = 1

        freedofs = np.arange(ndofs)
        freedofs = np.delete(freedofs, bcdofs)

        return Kstiff, Mmass, fint, freedofs, bcdofs

    def RayleighDamping(self, Kstiff, Mmass, zetas=[0.01, 0.01], OmegasS=[1, 2]):
        eigvals, eigvecs = eigs(Kstiff, k=10, M=Mmass, sigma=0)
        Omegas = np.sqrt(eigvals.real)
        Omega2 = Omegas[OmegasS[1]]
        Omega1 = Omegas[OmegasS[0]]

        if len(zetas) == 1:
            alpha, beta = 2 * Omega1 * zetas[0], 0
        elif len(zetas) == 2:
            beta = (2 * Omega2 * zetas[1] - 2 * Omega1 * zetas[0]) / (
                Omega2 ** 2 - Omega1 ** 2
            )
            alpha = 2 * zetas[0] * Omega1 - beta * Omega1 ** 2

        # Crayleigh = alpha * Mmass + beta * Kstiff

        return alpha, beta

    def assembleBeamElement(self, beam_mat_props, beam_cross, beam_Coords):
        Eyoung, nee, rho = beam_mat_props #[0]
        L = ln.norm(beam_Coords[0, :] - beam_Coords[1, :])
        MatrixT = self.beam_transform_matrix(
            beam_Coords[0, :], beam_Coords[1, :])

        A, I1, I2, I3 = beam_cross[0]
        f1, f2 = np.zeros(shape=(3)), np.zeros(shape=(3))
        elem_matrices = self.beam_mass_stiffness(
            L, A, I1, I2, I3, Eyoung, nee, rho, f1, f2
        )
        M, K, f = elem_matrices[0], elem_matrices[1], elem_matrices[2]
        Ke = (MatrixT.T.dot(K)).dot(MatrixT)
        Me = (MatrixT.T.dot(M)).dot(MatrixT)
        fe = MatrixT.T.dot(f)

        return Ke, Me, fe

    def elementglobaldofs(self, nodes):
        dofs = np.zeros(shape=(len(nodes) * 6))
        dofs[0: len(nodes) * 6: 6] = (nodes - 1) * 6
        dofs[1: len(nodes) * 6: 6] = (nodes - 1) * 6 + 1
        dofs[2: len(nodes) * 6: 6] = (nodes - 1) * 6 + 2
        dofs[3: len(nodes) * 6: 6] = (nodes - 1) * 6 + 3
        dofs[4: len(nodes) * 6: 6] = (nodes - 1) * 6 + 4
        dofs[5: len(nodes) * 6: 6] = (nodes - 1) * 6 + 5
        return dofs

    
    def beam_transform_matrix(self, node1, node2):
        ex = np.array([1, 0, 0])
        ey = np.array([0, 1, 0])
        ez = np.array([0, 0, 1])

        epsilons = self.beam_local_vectors(node1, node2)
        e1, e2, e3 = epsilons[0], epsilons[1], epsilons[2]

        Ti = np.array(
            [
                [np.dot(e1, ex), np.dot(e1, ey), np.dot(e1, ez)],
                [np.dot(e2, ex), np.dot(e2, ey), np.dot(e2, ez)],
                [np.dot(e3, ex), np.dot(e3, ey), np.dot(e3, ez)],
            ]
        )

        MatrixT = np.zeros(shape=(12, 12))
        MatrixT[0:3, 0:3] = Ti
        MatrixT[3:6, 3:6] = Ti
        MatrixT[6:9, 6:9] = Ti
        MatrixT[9:12, 9:12] = Ti

        # MatrixT = np.array([
        #    [Ti, np.zeros(shape=(3,3)), np.zeros(shape=(3,3)), np.zeros(shape=(3,3))],
        #    [np.zeros(shape=(3,3)), Ti, np.zeros(shape=(3,3)), np.zeros(shape=(3,3))],
        #    [np.zeros(shape=(3,3)), np.zeros(shape=(3,3)), Ti, np.zeros(shape=(3,3))],
        #    [np.zeros(shape=(3,3)), np.zeros(shape=(3,3)), Ti]
        # ])

        return MatrixT

    def AssembleSandZ(self, Mass):
        Zacceler = self.Excitation.SynthesizedAccelerogram
        ndofs = self.n_dofs
        originalnodes = self.Nodes.shape[0] - self.nl_link_elements.shape[0]
        angle = self.Excitation.angle
        loaddofsx = (np.arange(originalnodes)) * 6
        loaddofsy = (np.arange(originalnodes)) * 6 + 1

        alldofsx = (np.arange(ndofs))[0::6]
        alldofsy = (np.arange(ndofs))[1::6]
        alldofsz = (np.arange(ndofs))[2::6]

        Loaddofs = np.sort(np.hstack((loaddofsx, loaddofsy)))
        Alldofs = np.sort(np.hstack((alldofsx, alldofsy, alldofsz)))
        Mlumped = self.GetLumpedMass(Mass, Alldofs)
        Mlumped = np.diag(Mlumped)
        S = np.zeros(shape=(ndofs, 1))
        S[Loaddofs, 0] = Mlumped[Loaddofs]
        indexes = Loaddofs[::2]
        S[indexes, 0] = S[indexes, 0] * np.cos(angle)
        indexes = Loaddofs[1::2]
        S[indexes, 0] = S[indexes, 0] * np.sin(angle)

        return Zacceler, S

    def GetLumpedMass(self, Mass, translational_dofs):
        diagonal_coeff = np.diag(Mass)
        total_mass = np.sum(Mass)
        coeff_s = np.sum(diagonal_coeff[translational_dofs])
        lumped_diag = diagonal_coeff * total_mass / coeff_s
        Mlumped = np.diag(lumped_diag)
        return Mlumped

    def beam_local_vectors(self,node1, node2):
        ex = np.array([1, 0, 0])
        # ey=np.array([0,1,0])
        ez = np.array([0, 0, 1])

        e1 = (node2 - node1) / ln.norm(node2 - node1)
        is_all_zero = np.all((np.cross(e1, ez) == 0))
        if is_all_zero:
            e2 = np.cross(e1, ex) / ln.norm(np.cross(e1, ex))
        else:
            e2 = np.cross(e1, -ez) / ln.norm(np.cross(e1, -ez))

        e3 = np.cross(e1, e2)
        return [e1, e2, e3]

    def beam_mass_stiffness(self, L, A, I1, I2, I3, E, nee, rho, f1, f2):
        G = E / (2 * (1 + nee))

        Ka = np.array(
            [
                [E * A / L, 0, 0, 0, 0, 0],
                [0, 12 * E * I3 / (L ** 3), 0, 0, 0, 6 * E * I3 / (L ** 2)],
                [0, 0, 12 * E * I2 / (L ** 3), 0, -6 * E * I2 / (L ** 2), 0],
                [0, 0, 0, G * I1 / L, 0, 0],
                [0, 0, -6 * E * I2 / (L ** 2), 0, 4 * E * I2 / L, 0],
                [0, 6 * E * I3 / (L ** 2), 0, 0, 0, 4 * E * I3 / L],
            ]
        )

        Kb = np.array(
            [
                [-E * A / L, 0, 0, 0, 0, 0],
                [0, -12 * E * I3 / (L ** 3), 0, 0, 0, 6 * E * I3 / (L ** 2)],
                [0, 0, -12 * E * I2 / (L ** 3), 0, -6 * E * I2 / (L ** 2), 0],
                [0, 0, 0, -G * I1 / L, 0, 0],
                [0, 0, 6 * E * I2 / (L ** 2), 0, 2 * E * I2 / L, 0],
                [0, -6 * E * I3 / (L ** 2), 0, 0, 0, 2 * E * I3 / L],
            ]
        )

        Kc = np.array(
            [
                [E * A / L, 0, 0, 0, 0, 0],
                [0, 12 * E * I3 / (L ** 3), 0, 0, 0, -6 * E * I3 / (L ** 2)],
                [0, 0, 12 * E * I2 / (L ** 3), 0, 6 * E * I2 / (L ** 2), 0],
                [0, 0, 0, G * I1 / L, 0, 0],
                [0, 0, 6 * E * I2 / (L ** 2), 0, 4 * E * I2 / L, 0],
                [0, -6 * E * I3 / (L ** 2), 0, 0, 0, 4 * E * I3 / L],
            ]
        )

        K1 = np.zeros(shape=(12, 12))
        K1[0:6, 0:6] = Ka
        K1[0:6, 6:] = Kb
        K1[6:, 0:6] = Kb.T
        K1[6:, 6:] = Kc

        Ma = np.array(
            [
                [1 / 3, 0, 0, 0, 0, 0],
                [0, 13 / 35, 0, 0, 0, 11 * L / 210],
                [0, 0, 13 / 35, 0, -11 * L / 210, 0],
                [0, 0, 0, I1 / (3 * A), 0, 0],
                [0, 0, -11 * L / 210, 0, L ** 2 / 105, 0],
                [0, 11 * L / 210, 0, 0, 0, L ** 2 / 105],
            ]
        )

        Ma = Ma * rho * A * L

        Mb = np.array(
            [
                [1 / 6, 0, 0, 0, 0, 0],
                [0, 9 / 70, 0, 0, 0, -13 * L / 420],
                [0, 0, 9 / 70, 0, 13 * L / 420, 0],
                [0, 0, 0, -I1 / (6 * A), 0, 0],
                [0, 0, -13 * L / 420, 0, -(L ** 2) / 140, 0],
                [0, 13 * L / 420, 0, 0, 0, -(L ** 2) / 140],
            ]
        )

        Mb = Mb * rho * A * L

        Mc = np.array(
            [
                [1 / 3, 0, 0, 0, 0, 0],
                [0, 13 / 35, 0, 0, 0, -11 * L / 210],
                [0, 0, 13 / 35, 0, 11 * L / 210, 0],
                [0, 0, 0, I1 / (3 * A), 0, 0],
                [0, 0, 11 * L / 210, 0, L ** 2 / 105, 0],
                [0, -11 * L / 210, 0, 0, 0, L ** 2 / 105],
            ]
        )

        Mc = Mc * rho * A * L

        Ml = np.zeros(shape=(12, 12))
        Ml[0:6, 0:6] = Ma
        Ml[0:6, 6:] = Mb
        Ml[6:, 0:6] = Mb.T
        Ml[6:, 6:] = Mc

        fl = np.array(
            [
                [1 / 6 * L * (2 * f1[0] + f2[0])],
                [1 / 70 * L * (26 * f1[1] + 9 * f2[1])],
                [1 / 70 * L * (26 * f1[2] + 9 * f2[2])],
                [0],
                [-1 / 420 * (L ** 2) * (22 * f1[2] + 13 * f2[2])],
                [1 / 420 * (L ** 2) * (22 * f1[1] + 13 * f2[1])],
                [1 / 6 * L * (f1[0] + 2 * f2[0])],
                [1 / 70 * L * (9 * f1[1] + 26 * f2[1])],
                [1 / 70 * L * (9 * f1[2] + 26 * f2[2])],
                [0],
                [1 / 420 * (L ** 2) * (13 * f1[2] + 22 * f2[2])],
                [-1 / 420 * (L ** 2) * (13 * f1[1] + 22 * f2[1])],
            ]
        )

        return [Ml, K1, fl]
