ReadMe file for orientation purposes 

The ExampleMainFile notebook simulates the response of the two-story frame benchmark with Bouc-Wen hysteretic 
links under ground motion excitation.

Each section of the notebook described the definition of a key component on assembling the two-story frame.
Short descriptions and orientation remarks are also provided inside the notebook.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Input FIle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The input file class assembles the necessary arrays to define the geometrical configuration of the 
two-story frame. It also offers the possibility of assembling nonlinear hysteretic links on the frame,
modeled following the Bouc-Wen hysteretic formulation. The more advanced possibilities for spring and plate
elements are not yet available and can be found only in the MATLAB version of the benchmark.

The example input file provided in the notebook follows the InputFileLinksAll file of the MATLAB version, 
where hysteretic links are located in all horizontal beams and the basement and first-story columns of the two-story frame.
In the first story columns there is a link at the bottom and one at the top.

The main attributes of the class object are explained here:

.Nodes
Node coordinates of the frame
Nodes with a nonlinear hysteretic link are inserted here twice for modeling purposes. 
For example, if a Bouc Wen link exists between node 1 and 2, these two nodes must be defined twice as
different entities.

.BeamElements
Connectivity matrix for beam elements
Only linear beam elements are assembled here. The notation is explained on the file.

.nl_link_elements
Connectivity matrix for nonlinear links
Nonlinear links are assembled in this matrix using the virtual duplicate nodes assembled in MODEL.nodes
For example, if two links exist between node 1 and 2 then:
-We need to define node 11 as a duplicate of node 1
-We need to define node 12 as a duplicate of node 2
-We assemble the beam elements 11-12
-We assemble the nonlinear links 1-11 and 12-2 etc.

.nl_link_flags
Dofs where the nonlinear links are activated
The notation is explained on the notebook.

.NodalDisplacements
Array of nodal displacements. Boundary conditions are defined here. 
Each row corresponds to a node. The first six columns are flags to indicate whether a constraint is applied in the corresponding dof.
The last six columns are the values of the applied displacements. The exact notation is also described in the notebook. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Material Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Material class:
A class containing the material parameters of the model and the respective
cross sections.
The attributes are self-explainable.
Eyoung denotes the Young modulus
nee denotes the poisson ratio
rho denotes the density of the material.
In addition, the cross section properties are hardcoded inside the class.
Eg cross section area A, moments of inertia etc.
Currently only one cross section is provided as reference.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hysteretic Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Regarding the Bouc-Wen model, all links assume the same parametric model.
The advanced option to modify certain links already available in the MATLAB version of the benchmark
is not provided here. The user is able to follow the MATLAB notation as a guideline and adjust the notebook
if needed. 

As already described, a Bouc-Wen model is injected at the links indicated at
InputFile.nl_link_elements and into the dofs indicated at the InputFile.nl_link_flags matrix.

The BoucWenModel class controls the parameters of the hysteretic links. 
It contains the parameters for the hysteretic mapping and a method to evaluate
the forcing on all links given the displacements, the history variables and a
degree of freedom selection matrix to properly assembly the incremental displacements
for each link. The updated version includes stiffness degradation and strength deterioration
controlled by parameters deltav and delta n respectively.

The physics of the Bouc-Wen model are explained on the attached pdf description file.
As mentioned, on the current template all links assume the same  displacement-hysteretic force mapping.

Parameters simplified explanation: 
Amplitude of hysteresis is controled by parameter AmpBW and Alpha
Shape of hysteretic curve by parameters beta and gamma.
No normalization is assumed.
Parameters deltav and deltan control stiffness degradation and strength deterioration.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Excitation and Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The SimulationParameters class defines the ground motion excitation signal.

Several alternatives are provided, based on the "type" input variable:
"ExampleA" / "ExampleB" are pre-defined example signals 
"multisine" constructs a multisine signal based on the user defined spectral and temporal components
"whitenoise" constructs a white-noise like signal based on the user defined spectral and temporal components

To construct the excitation signal the following inputs are needed:
fmin: minimum frequency component
fmax: maximum frequency component
fs: sampling frequency
AmpF: amplitude coefficient
Ncomponents: Number of frequency components (multisine) or Number of samples (white noise)
fstep: Variable to take every fstep frequencies to construct the multisine signal (fmin:fstep:fmax),
if None the step value is defined using the Ncomponents. If fstep=0 and Ncomponents=None 10 random
frequency components are used to construct the signal.

The angle of motion is defined by the input variable angle. 
Currently only an in-plane(x-y) ground motion is provided.
The user can follow the MATLAB version's notation if a more complex motion is desired.
Thus, the angle is formulated with respect to the x longitudinal axis. This implies
that the excitation on the x dofs will be multiplied by cos(Angle) whereas
for the y direction by sin(Angle). The nodal excitation is assembled by 
multiplying the accelerogram with the (diagonal) nodal mass contributions.

The user can also provide an input excitation by setting the *type* variable to "user" and 
providing the signal as an numpy 1D array in the *UserSignal* entry.

The initial conditions are defined in the Main File. If not, zero initial conditions are assumed.
The notations is ICs["Displacements"]=, ICs["Velocities"]= and require a vector with ndof entries.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Assembly definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

As already explained:
InputFile takes care of the input file and the geometrical configuration definition
Material defines material properties
BoucWenModel assembles the hysteretic model of the links
SimulationParameters contains the excitation signal and the initial conditions.

Assembly takes all previous components as inputs and assembles the model to be
simulated. This implies that the class assembles the initial matrices and vectors
based on the respective equation of motion (K,M,C,F,R) and defines all relevant
initializations or parameters for the simulation.

Important: DAMPING!
Damping is defined based on the Rayleigh assumption.
You can directly give the coefficients or define the zeta rations to perform
modal analysis and extract the coefficients. There is a dedicated method for that
inside the class.The user can define manual damping by setting 
SimulationParameters.ICs["zeta"]=None and SimulationParameters.ICs["OmegaIndexes"]=[alpha, beta]

ATTENTION:
Rayleigh_b multiplies the stiffness matrix and Rayleigh_a the mass matrix.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Post processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The degree of freedom notation is:
1 - x - Axial direction
2 - y - shear in plane direction
3 - z - height - shear out of plane
4 - Torsion
5(I2),6(I3) - Rotational degrees of freedom 

Axis x is defined as the longitudinal/axial axis and assumes
the direction of the largest span of the frame.

This implies that every node has 6 corresponding rows in the results matrices.
Eg. In U or HistU the 1st row is the x dof response of the first node.
For the local system, you can inspect the notation by checking the K contributions in the beam_local_stiffness_mass_0
function. The notation follows an equivalent formulation with rotation applied when needed.


Format of the Results object:
Results["Disps"] : Displacement Time histories
Results["Velocities"] : Velocities Time histories
Results["AccelerationsA"]  : Acceleration Time histories
Every row corresponds to a degree of freedom (spatial discretization) and every column to a timestep (time integration)

Results["HystereticParamsU"] = Displacement time history of nonlinear links.
This is not the hysteretic variable zeta.
Results["HystereticParamsR"] = Hysteretic forcing time history of nonlinear links.
