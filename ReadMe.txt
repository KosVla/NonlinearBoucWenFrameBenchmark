ReadMe file for orientation purposes 

The MainScript* simulates the response of a frame under a dynamic excitation.
Example Main Scripts:
ExampleMainFile_AllPossibleInputs
ExampleMainFile_Clean (clean version of the same input file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INPUT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The input file offers the possibility to assemble nonlinear hysteretic links on the frame, 
modeled after the Bouc Wen hysteresis formulation. There are also advanced possibilities for
spring elements, plate elements etc. Going through the InputFile you will find that
each matrix/section is explained in comments.
There is a dedicated folder for input files:
InputFileLinksFirst: BW links in all beam connections and the columns of the basement.
In the columns the link is located at the top.
InputFileLinksFirst_DifferentBW: BW links in all beam connections and the columns of the basement
				and additionally demonstrates how to change the BW parameters of individual links
InputFileLinksBeams: BW links only in beam connections.
InputFileLinksAll: Hysteretic links are located in all horizontal 
beams and columns of the two-story frame. In the columns there is a link 
at the bottom and one at the top of the column.

The input file creates a MODEL struct with necessary parameters for the frame simulation.
*The hysteretic parameter allocation is handled in the GetHystereticParams function!!!

Going through the components of the MODEL struct:

MODEL.nodes : Node coordinates of the frame
Nodes with a nonlinear hysteretic link are inserted here twice for modeling purposes. 
For example, if a Bouc Wen link exists between node 1 and 2, these two nodes must be defined twice as
different entities.

MODEL.beam_elements : Connectivity matrix for beam elements
Only linear beam elements are assembled here. The notation is explained on the file.

MODEL.nl_link_elements : Connectivity matrix for nonlinear links
Nonlinear links are assembled in this matrix using the virtual duplicate nodes assembled in MODEL.nodes
For example, if two links exist between node 1 and 2 then:
-We need to define node 11 as a duplicate of node 1
-We need to define node 12 as a duplicate of node 2
-We assemble the beam elements 11-12
-We assemble the nonlinear links 1-11 and 12-2 etc.

MODEL.material_properties : Material properties
The notation is explained on the file.

MODEL.beam_material_properties: Material properties of beam elements
The notation is explained on the file.

MODEL.nl_links_alternate:
Nonlinear links with different BW parameters than the Input struct defined on the Main file
The notation is explained on the file.

MODEL.nl_link_flags : Dofs where the nonlinear links are activated
The notation is explained on the file.
By setting a flag equal to zero the respective Bouc-Wen link will behave linearly

MODEL.cross_sections : Cross section properties
MODEL.beam_cross_sections: Selection vector of the beam cross sections
The notation is explained on the file.

MODEL.beam_loads : Loads on beam elements
The notation is explained on the file.

MODEL.nodal_loads : Load on nodes
The notation is explained on the file.

MODEL.nodal_displacements : BCs and nodal displacements
The notation is explained on the file.

MODEL.springs : Additional spring elements
The notation is explained on the file.

MODEL.masses : Additional lumped masses on nodes
The notation is explained on the file.

MODEL.dyn.dt : time step
MODEL.dyn.nt : duration of the analysis
MODEL.dyn.a : Rayleigh damping coefficient alpha (Mass proportional)
MODEL.dyn.b : Rayleigh damping coefficient beta (Stiffness proportional)
*There is also the option of calculating those coefficients by giving the 
zeta ratios and the Omegas desired in the main file. If you don't give any zeta
then alpha and beta are defined in the input file. 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
The AutomaticMesh function generates the nodal coordinates and connectivity,
along with the respective matrices for properties definition, for a user
specified number of floors, frames in x, frames in y and the respective dimensions.
This allows the user to create a frame based on their needs.
The documentation of the function describes the input parameters in detail. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Excitation and Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The excitation and/or the initial conditions are defined in the Main File.
Two examples are provided for the excitation, along with a multisine constructor
for parametrized excitation signals. The notation is described in the description
file provided as well. On the clean version there is a dedicated function for handling
the assembly of the excitation vector. Parameters that control the excitation
in general are the following:

Input.Angle
Multiple alternatives are provided. 
The user an have a look on the ExampleMainFile_AllPossibleInputs
for on averview.
This variable stores the angle of excitation.
If it is a single value, it implies a xy plane excitation
with a 2D angle in the xy plane with axis x, as defined in the description.
If it is a 3D angle, the user should follow the notation below:
[ (Angle of the 3D vector with respect to z axis) (Angle of the XY plane projection with axis x) Quadrant of 3D vector]
Finally, the user can also define the coordinates/weights of the 3D excitation vector by providing
[Coordx Coordy Coordz 0]. The software then assumes that the ground motion is multiplied by Coordx in the x axis and so on.
For example:
2D angle in xy plane with axis x: Angle=pi/4;
3D angle with the following notation:
[Angle(3D vector with z axis), Angle(XY plane projection with x), Quadrant]
eg. Angle=[ pi/2 pi/4 1] is equivalent to the first alternative.
3D vector [Coordx Coordy Coordz 0]
Angle=[cos(pi/4) sin(pi/4) 0 0] is equivalent to the first alternative.


The excitation signal is assumed to be a ground motion acceleration.
It is stored in Input.SynthesizedAccelerogram.
So, if the user desires to define a specific accelerogram they can change this entry.
Nodal excitation can be defined from the input file only!!!
Regarding the parametrization of the excitation, the following notation is followed:
P : number of excitation periods.
N : number of points per period.
upsamp: upsampling factor
fmin / fmax: excitation bandwidth.
AmpF or A: excitation amplitude coefficient.

The nodal excitation is assembled by multiplying the accelerogram with the nodal mass contributions from the lumped mass matrix.
However, the integration uses the consistent mass matrix!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hysteretic Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Regarding the Bouc-Wen model, all links assume the same parametric model.
The ones with different parameters should be defined in the Input file under 
MODEL.nl_links_alternate following the respective notation explained previously.

As already described, a Bouc-Wen model is injected at the links indicated at
MODEL.nl_link_elements and into the dofs indicated at the MODEL.nl_link_flags matrix.
These matrices can be adjusted in the input file to refer to group of or individual links as well.

The physics of the Bouc-Wen model are explained on the attached description.
The parameters of the Bouc-Wen are defined in the main file and are assembled in the
GetHystereticParams function inside the BoucWenRun function. In this function we also
initialize our initial conditions matrices, load matrices and our history struct.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
The Bouc-Wen parameters of the links are assembled in the GetDofsLoop function.
The output of this function is a matrix that contains the starting and end dofs of each hysteretic link
in the first two columns and then the rest 9 columns contain the individual BW parameters for each link and 
each degree of freedom of the link.
If the degree of freedom is not activated based on the nl_link_flags matrix, the respective parameters of the matrix
are modified to enforce a linear behavior.
If the element is contained on the nl_links_alternate input matrix then the parameters of the matrix are modifies 
to correspond to the input ones (different than the other elements).

The Bouc-Wen model evaluation happens at bw_Ndall.
This function takes as input vectors with the parameters of all hysteretic links assembled properly, based on the
GetDofsLoop function. Thus, the function is vectorized and evaluates all links at once. 
The user can modify this function as desired, as long as the vectorized nature of the input and outputs is retained.

Along with the vectorized version of the code, a dedicated core folder (coreBW_SDOF) is provided that 
assembles each degree of freedom of the hysteretic links separately with a SDOF Bouc-Wen function (bw_1d). 
This allows the user to distinguish between links and implement different models if desired.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Notation of MODEL.DofsLoop:
Each row refers to a degree of freedom (activated or not)
[StartDof EndDof a k A b g w AmpBW deltav deltan] 


Bouc-Wen hysteresis model parameters - Notation follows description pdf
bw_a=a / Alpha=A / N=w / Beta=b / Gamma=g / deltav=dv / deltan=dnu / bw_k=k
AmpBW additional amplitude parameter for BW. Multiplies only hysteretic term to magnify influence
Simplified explanation:
Aplitude of hysteresis is controled by parameter AmpBW, Alpha
Shape and smoothness of hysteretic curve by parameters beta and gamma.
Parameters deltav and deltan control stiffness degradation and strength deterioration.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Post processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results struct : MODEL

Based on Figure 1 in the description, the degree of freedom notation for the global results matrices is:
1 - x 
2 - y 
3 - z 
4 - Torsion
5(I2),6(I3) - Rotational degrees of freedom 

This implies that every node has 6 corresponding rows in the results matrices.
Eg. In U or HistU the 1st row is the x dof response of the first node.
For the local system, you can inspect the notation by checking the K contributions in the beam_local_stiffness_mass_0
function. The notation follows an equivalent formulation with rotation applied when needed.

MODEL.Hist* : Time histories of nonlinear links.
MODEL.HistU = Displacement time history of nonlinear links
MODEL.HistR = Hysteretic forcing time history ofnonlinear links
Each row corresponds to a degree of freedom
E.g. to find the first dof of link j you need to evaluate the following:
(number of entries equal to 1 in MODEL.nl_link_flags)*(j-1)+1

MODEL.U : Displacement Time histories
MODEL.V : Velocities Time histories
MODEL.A : Acceleration Time histories

MODEL.K : Stiffness matrix
MODEL.M : Mass matrix
MODEL.f : External excitation at last timestep
MODEL.fint : Internal forces at last timestep
MODEL.Rmatrix : External excitation time history
MODEL.u : Displacements at last timestep

Function plot_model :
Visualizes the deformed configuration of the model in 3D. By changing MODEL.u you can visualize
a different timestep (not the last one) or any other meaningful quantity.