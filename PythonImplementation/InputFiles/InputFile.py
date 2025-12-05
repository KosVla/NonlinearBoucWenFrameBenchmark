# -*- coding: utf-8 -*-

__author__ = "Konstantinos Vlachas"
__email__ = "vlachask@ibk.baug.ethz.ch"


import numpy as np

def AssembleFrame(MODEL, InputFile='FloorsExample', lengthx=7.50, lengthy=5.00, heightz=3.2):
    match InputFile:
        case 'FloorsExample':
            AssembleFrame_FloorsExample(MODEL, lengthx, lengthy, heightz)
        case 'BeamsExample':
            AssembleFrame_BeamsExample(MODEL)   

def AssembleFrame_BeamsExample(MODEL, lengthx=7.50, lengthy=5.00, heightz=3.2):
    lx,ly, lz = lengthx, lengthy, heightz
    
    MODEL.Nodes =np.array([
        [0,  0,  0], #0
        [lx, 0,  0], #1
        [2*lx, 0, 0], #2
        [2*lx, ly, 0], #3
        [lx, ly, 0], #4
        [0,  ly, 0], #5
        [0,  0,  lz], #6
        [lx, 0,  lz], #7
        [2*lx, 0, lz], #8
        [2*lx, ly, lz], #9
        [lx, ly, lz], #10
        [0,  ly, lz], #11
        [0,  0, 2*lz], #12
        [lx, 0,  2*lz], #13
        [2*lx, 0, 2*lz], #14
        [2*lx, ly, 2*lz], #15
        [lx, ly, 2*lz], #16
        [0,  ly, 2*lz], #17
        [0,  0,  0], #18
        [lx, 0,  0], #19
        [2*lx, 0, 0], #20
        [2*lx, ly, 0], #21
        [lx, ly, 0], #22
        [0,  ly, 0], #23
        [0,  0,  lz],#24
        [0,  0,  lz],#25
        [lx, 0,  lz],#26
        [lx, 0,  lz],#27
        [lx, 0,  lz],#29
        [2*lx, 0, lz],#31
        [2*lx, 0, lz],#33
        [2*lx, ly, lz], #34
        [2*lx, ly, lz], #36
        [lx, ly, lz], #37
        [lx, ly, lz], #38
        [lx, ly, lz], #40
        [0,  ly, lz],#41
        [0,  ly, lz],#43
        [0,  0,  2*lz], #44
        [0,  0,  2*lz], #45
        [lx, 0,  2*lz], #46
        [lx, 0,  2*lz], #47
        [lx, 0,  2*lz], #48
        [2*lx, 0,  2*lz], #49
        [2*lx, 0,  2*lz], #50
        [2*lx, ly, 2*lz], #51
        [2*lx, ly, 2*lz], #52
        [lx, ly, 2*lz], #53
        [lx, ly, 2*lz], #54
        [lx, ly, 2*lz], #55
        [0,  ly, 2*lz], #56
        [0,  ly, 2*lz], #57         
    ])
            
    MODEL.BeamElements = np.array([
        [0, 1, 7],
        [0, 2, 8],
        [0, 3, 9],
        [0, 4, 10],
        [0, 5, 11],
        [0, 6, 12],
        [0, 7, 13],
        [0, 8, 14],
        [0, 9, 15],
        [0, 10, 16],
        [0, 11, 17],
        [0, 12, 18],
        [0, 38, 25],
        [0, 26, 27],
        [0, 29, 30],
        [0, 28, 35],
        [0, 31, 32],
        [0, 33, 34],
        [0, 36, 37],
        [0, 52, 39],
        [0, 40, 41],
        [0, 43, 44],
        [0, 42, 49],
        [0, 45, 46],
        [0, 47, 48],
        [0, 50, 51],
    ])

    nl_link_elements = np.array([
        [19, 1, 0],
        [20, 2, 0],
        [21, 3, 0],
        [22, 4, 0],
        [23, 5, 0],
        [24, 6, 0],
        #Beams first floor
        [25, 7, 2],
        [7 ,26, 2],
        [27, 8, 2],
        [8, 28, 2],
        [8, 29, 2],                
        [30, 9, 2],
        [9, 31, 2],
        [32, 10, 2],
        [10, 33, 2],                
        [34, 11, 2],                
        [35, 11, 2],                
        [11, 36, 2],
        [37, 12, 2],
        [12, 38, 2],
        #Beams second floor
        [39, 13, 3],
        [13, 40, 3],
        [41, 14, 3],
        [14, 43, 3],
        [14, 42, 3],
        [44, 15, 3],
        [15, 45, 3],
        [46, 16, 3],
        [16, 47, 3],
        [48, 17, 3],
        [49, 17, 3],
        [17, 50, 3],
        [51, 18, 3],
        [18, 52, 3],
    ])         
            
    MODEL.nl_link_elements = nl_link_elements

    MODEL.nl_link_flags = np.ones(shape=(nl_link_elements.shape[0], 6))

    ground = np.array([19, 20, 21, 22, 23, 24])

    NodalDisplacements = np.zeros(shape=(ground.shape[0], 13))
    NodalDisplacements[:, 0] = ground
    NodalDisplacements[:, 1:7] = 1
    MODEL.NodalDisplacements = NodalDisplacements

    # Define Initial Conditions
    ndim = np.shape(MODEL.Nodes)[0] * 6
    ICs = dict()
    ICs["Disps"] = np.zeros((ndim, 1))
    ICs["Velocities"] = np.zeros((ndim, 1))
    ICs["Accelerations"] = np.zeros((ndim, 1))
    ICs["zeta"] = None
    #ICs["OmegaIndexes"] = [0, 1]
    MODEL.ICs = ICs
    MODEL.n_dofs = ndim

def AssembleFrame_FloorsExample(MODEL, lengthx=7.50, lengthy=5.00, heightz=3.2):
    lx,ly, lz = lengthx, lengthy, heightz
    
    MODEL.Nodes =np.array([
        [0,  0,  0], #0
        [lx, 0,  0], #1
        [2*lx, 0, 0], #2
        [2*lx, ly, 0], #3
        [lx, ly, 0], #4
        [0,  ly, 0], #5
        [0,  0,  lz], #6
        [lx, 0,  lz], #7
        [2*lx, 0, lz], #8
        [2*lx, ly, lz], #9
        [lx, ly, lz], #10
        [0,  ly, lz], #11
        [0,  0, 2*lz], #12
        [lx, 0,  2*lz], #13
        [2*lx, 0, 2*lz], #14
        [2*lx, ly, 2*lz], #15
        [lx, ly, 2*lz], #16
        [0,  ly, 2*lz], #17
        [0,  0,  0], #18
        [lx, 0,  0], #19
        [2*lx, 0, 0], #20
        [2*lx, ly, 0], #21
        [lx, ly, 0], #22
        [0,  ly, 0], #23
        [0,  0,  lz],#24
        [0,  0,  lz],#25
        [0,  0,  lz],#26
        [lx, 0,  lz],#27
        [lx, 0,  lz],#28
        [lx, 0,  lz],#29
        [lx, 0,  lz],#30
        [2*lx, 0, lz],#31
        [2*lx, 0, lz],#32
        [2*lx, 0, lz],#33
        [2*lx, ly, lz], #34
        [2*lx, ly, lz], #35
        [2*lx, ly, lz], #36
        [lx, ly, lz], #37
        [lx, ly, lz], #38
        [lx, ly, lz], #39
        [lx, ly, lz], #40
        [0,  ly, lz],#41
        [0,  ly, lz],#42
        [0,  ly, lz],#43
        [0,  0,  2*lz], #44
        [0,  0,  2*lz], #45
        [lx, 0,  2*lz], #46
        [lx, 0,  2*lz], #47
        [lx, 0,  2*lz], #48
        [2*lx, 0,  2*lz], #49
        [2*lx, 0,  2*lz], #50
        [2*lx, ly, 2*lz], #51
        [2*lx, ly, 2*lz], #52
        [lx, ly, 2*lz], #53
        [lx, ly, 2*lz], #54
        [lx, ly, 2*lz], #55
        [0,  ly, 2*lz], #56
        [0,  ly, 2*lz], #57         
    ])
            
    MODEL.BeamElements = np.array([
        #Columns basement
        [0, 1, 27],
        [0, 2, 31],
        [0, 3, 34],
        [0, 4, 37],
        [0, 5, 40],
        [0, 6, 44],
        #Beams first floor
        [1, 26, 28],
        [1, 30, 32],
        [1, 33, 35],
        [1, 36, 41],
        [1, 39, 42],
        [1, 43, 25],
        [1, 29, 38],
        #Columns first floor
        [2, 7, 13],
        [2, 8, 14],
        [2, 9, 15],
        [2, 10, 16],
        [2, 11, 17],
        [2, 12, 18],
        #Beams second floor
        [3, 46, 47],
        [3, 49, 50],
        [3, 51, 52],
        [3, 53, 54],
        [3, 56, 57],
        [3, 58, 45],
        [3, 48, 55],
    ])

    nl_link_elements = np.array([
        #Links on the basement nodes
        [19, 1, 0],
        [20, 2, 0],
        [21, 3, 0],
        [22, 4, 0],
        [23, 5, 0],
        [24, 6, 0],
        #Columns first floor
        [27, 7, 1],
        [31, 8, 1],
        [34, 9, 1],
        [37, 10, 1],
        [40, 11, 1],
        [44, 12, 1],
        #Beams first floor
        [25, 7, 2],
        [7 ,26, 2],
        [28, 8, 2],
        [29, 8, 2],
        [8, 30, 2],                
        [32, 9, 2],
        [9, 33, 2],
        [35, 10, 2],
        [10, 36, 2],                
        [38, 11, 2],                
        [11, 39, 2],                
        [41, 11, 2],
        [42, 12, 2],
        [12, 43, 2],
        #Beams second floor
        [45, 13, 3],
        [13, 46, 3],
        [47, 14, 3],
        [14, 48, 3],
        [14, 49, 3],
        [50, 15, 3],
        [15, 51, 3],
        [52, 16, 3],
        [16, 53, 3],
        [54, 17, 3],
        [55, 17, 3],
        [17, 56, 3],
        [57, 18, 3],
        [18, 58, 3],
    ])         
            
    MODEL.nl_link_elements = nl_link_elements

    MODEL.nl_link_flags = np.ones(shape=(nl_link_elements.shape[0], 6))

    ground = np.array([19, 20, 21, 22, 23, 24])

    NodalDisplacements = np.zeros(shape=(ground.shape[0], 13))
    NodalDisplacements[:, 0] = ground
    NodalDisplacements[:, 1:7] = 1
    MODEL.NodalDisplacements = NodalDisplacements

    # Define Initial Conditions
    ndim = np.shape(MODEL.Nodes)[0] * 6
    ICs = dict()
    ICs["Disps"] = np.zeros((ndim, 1))
    ICs["Velocities"] = np.zeros((ndim, 1))
    ICs["Accelerations"] = np.zeros((ndim, 1))
    ICs["zeta"] = [0.02, 0.02]
    ICs["OmegaIndexes"] = [0, 1]
    MODEL.ICs = ICs
    MODEL.n_dofs = ndim
        
            
class Material:
    def __init__(self, Eyoung, nee=0.3, rho=8000):
        
        if isinstance(Eyoung, list):
            nee, rho = [nee for x in Eyoung], [rho for x in Eyoung]
            mat_props = np.array(list(zip(Eyoung, nee, rho)))
        else:
            EyoungL, nee, rho =[Eyoung for x in [1,1,1,1]], [nee for x in [1,1,1,1]], [rho for x in [1,1,1,1]]
            mat_props = np.array(list(zip(EyoungL, nee, rho)))

        self.mat_properties=mat_props

        # Steel HEA 200 cross section
        A, I1, I2, I3 = 5383e-06, 204.3e-09, 36.92e-06, 13.36e-06

        self.cross_sections = np.array([[A, I1, I2, I3]])