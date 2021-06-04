******************************************************************************************************************
ReadMe file for orientation purposes
Please refer to the description of the benchmark for more details on how the datasets were created
******************************************************************************************************************

****Configurations Setup****

*ConfigurationA:
Mechanical damage is induced in the healthy state of the frame.
The healthy state is represented by the linearized model of the frame. The
Bouc-Wen links thus produce a linear behavior using the following set-up:
a= 1.0, k= 6e7, AmpBW=0.0 (The rest of the parameters do not matter).

The induced mechanical damage is represented by the hysteretic Bouc-Wen 
links simulated with the following parametric set-up:
a= 0.5276, k= 5.4236e7, beta=12.6, gama=5.7, deltah=deltav=0

For each input signal both the healthy(MODELH) and the damaged (MODELD) 
model are simulated. Their displacement output is stored on the respective
workspace. Based on the relative error norm(MODELD.U-MODELH.U)/norm(MODELH.U)
each pair of healthy and damaged simulations is characterized as
hard/medium/easy (errors: <10/10-20%/>20%)

*ConfigurationB:
Damage growth is induced in the initial state of the frame representing
a deterioration scenario. Here, the initial/healthy state is represented
the frame with activated links and the following set-up:
beta= 3.54,gama= 6.39, a= 1.0, k= 6e7, AmpBW=1.0, deltah=deltav=0
The damage/deterioration scenario is represented by the hysteretic Bouc-Wen 
links simulated with the following parametric set-up:
beta= 3.54,gama= 6.39, a= 1.0, k= 6e7, AmpBW=1.0, (same as the initial)
deltav = 2.57, deltan=1.42.

For each input signal both the initial(MODELH) and the damaged (MODELD) 
model are simulated. Their displacement output is stored on the respective
workspace. Based on the relative error norm(MODELD.U-MODELH.U)/norm(MODELH.U)
each simulation is characterized as hard/medium/easy (errors: <10/10-20%/>20%)

*ConfigurationC:
Model calibration scenario during operation stages
Here, the initial/healthy state is represented by the frame with activated
links and the following set-up:
 beta= 4.67,gama= 1.21, a= 0.10, k= 1e8, AmpBW=1.0, deltah=deltav=0
The model calibration part requires estimating the parameters of the model
for a changed state, represented by a new hysteretic links set-up:
beta= 30.87,gama= -40.21, a= 0.10, k= 1e8, AmpBW=1.0, deltah=deltav=0

For each input signal both the healthy(MODELH) and the calibrated (MODELD) 
model are simulated. Their displacement output is stored on the respective
workspace. Based on the relative error norm(MODELD.U-MODELH.U)/norm(MODELH.U)
each simulation is characterized as hard/medium/easy (errors: <10/10-20%/>20%)

*ConfigurationD:
Localized damage concentrated on all links of the first floor (beams+columns)
Here, the initial/healthy state is represented by the frame with activated
links and the following set-up:
beta=21.82,gama=9.71, a= 0.6235, k= 1e8, AmpBW=1.0, deltah=deltav=0
The induced mechanical damage is represented by the hysteretic Bouc-Wen 
links (first floor only) simulated with the following parametric set-up:
beta=21.82,gama=9.71, a=0.1236, k=7.7234e7, AmpBW=1.0, deltah=deltav=0

For each input signal both the initial(MODELH) and the calibrated (MODELD) 
model are simulated. Their displacement output is stored on the respective
workspace. Based on the relative error norm(MODELD.U-MODELH.U)/norm(MODELH.U)
each simulation is characterized as hard/medium (errors: <10/10-20%)

*ConfigurationE:
Damage growth is induced in the initial state of the frame representing
a localized deterioration scenario (first floor links only).
Here, the initial/healthy state is represented the frame with activated 
links and the following set-up::
beta=16.32,gama= -16.32, a= 0.10, k= 5e7, AmpBW=1.0, deltah=deltav=0
The damage/deterioration scenario is represented by the hysteretic Bouc-Wen 
links simulated with the following parametric set-up (first floor links only):
beta=16.32,gama= -16.32, a= 0.10, k= 5e7, AmpBW=1.0, deltan=4.85,
deltav=2.31

For each input signal both the initial(MODELH) and the damaged (MODELD) 
model are simulated. Their displacement output is stored on the respective
workspace. Based on the relative error norm(MODELD.U-MODELH.U)/norm(MODELH.U)
each simulation is characterized as hard/medium (errors: <10/10-20%)


****Naming convention****

Example: Configuration(A)_(Easy)Task

*Configuration(A or B or C or D or E) 
Refers to the standardized cases described in the detailed report of the 
benchmark simulator and in this ReadMe file.
Five configurations are provided, each one produced with a different 
standardized parametric set.

*(Easy or Medium or Hard)Task
Refers to the intensity of the nonlinear behavior, or more precisely,
to how different the actual signal is from the reference configuration due 
to the different parametric setup. For details check the configuration description
above. For example:
-Configuration A: Refers to the intensity of the nonlinear effects 
with respect to the healthy state of the system (No hysteretic links, linear rigid connections).
Thus, in the Easy version the actual signal has high errors compared against
the healthy state and it is relatively easy to distinguish that something is off,
whereas the hard case has lower errors and is closer to the healthy response.
-Configuration BC and E: Refers to the initial state of the model, as presented in Table 3. These two cases
feature a shear frame with Bouc-Wen links as an initial model and the parametrization of the links changes for some reason,
representing damage or similar phenomena. Thus, the easy/medium/hard tasks refers to how different the actual signal is
compared against the initial configuration. More details in Table 3.

*(Test or Train or Validation)Case_No*
Each configuration features three tasks levels and each task features three datasets:
A training, a testing and a validation one. The user is free to shuffle those, if desired.
However, the Configuration used should be clearly reported for comparison and reference
purposes. 


****Workspace variables****
Every .zip file contains workspaces for the signals and for the 
parametrization used to create the responses.
The parametrization refers to the excitation characteristics, as described
in the report provided (minimum frequency
of multisine excitation, amplitude coefficient, angle of motion in xy plane etc).
Each one of the response workspaces contains the Results struct. 
There, the following fields are provided:
-OutputU: Output displacement response of the mdof oscillator (all degrees of freedom)
-OutputUH: Output displacement response of the mdof oscillator (initial/healthy state)
-Input: Struct containing the model parametrization
(see detailed report for notation, the input excitation is contained here
under Input.SynthesizedAccelerogram)
-A, fmin,fmax, Angle, etc: Excitation parameters.


