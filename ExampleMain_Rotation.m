%% Example main file to evaluate the response of two-story shear frame 
% with BW links 
%
% Modified Features compared to Version 0.0:
% -Bouc-Wen links are introduced only on the rotational degrees of freedom
% -Bending and shear degrees of freedom are now coupled
% -The option to assemble the Bouc-Wen bw_k coefficient from the EI coefficient
% of the respective beam element is provided
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.
%and
% J. Noel and M. Schoukens, 
% Hysteretic benchmark with a dynamic nonlinearity,
% in Workshop on non-linear system identification benchmarks, 2016, pp. 7–14

%% Clear workspace and load folders
clc;clear;
restoredefaultpath
addpath coreRotationOnly
addpath InputFiles

%% Define Input
% Definition of input files for the analysis

%Select Model 

% Two story frame, only rotational degrees of freedom are activated
% MODEL = InputFileLinksAll_RotationOnly();  

%AutomaticMeshColumns examples - Links only on the columns
%NoofFloors excludes the basement, so floors additional to the basement
% MODEL = AutomaticMesh(NoOfFloors,NoOfFramesx,NoOfFramesy,dimensions); 
MODEL = AutomaticMeshColumns(6,2,1,[7.50,5.0,3.2]);  

ndim=6*numel(MODEL.nodes(:,1));
%Plot undeformed state of the model and visualize beam orientation
% plot_model( MODEL, 0 )


%Define parametric dependencies

%Bouc-Wen hysteresis model parameters - Notation follows description pdf
%bw_a=a / Alpha=A / N=w / Beta=b / Gamma=g / deltav=dv / deltan=dnu
%bw_k=k
Input.Alpha=1.00; Input.N=1; Input.bw_a = 0.10;
Input.bw_k = 1e7;
Input.Beta=3; Input.Gamma=2; Input.deltav = 0; Input.deltan=0;
%Additional amplitude parameter for BW. Multiplies only hysteretic term to
%magnify influence
Input.AmpBW=1.0; 

%Option to assemble the Bouc-Wen bw_k coefficient from the EI coefficient
%of the respective beam element is provided (Input.BWkfromMAT=true)
Input.BWkfromMAT=true;

%Initial conditions (displacements/velocities/accelerations)
Input.u0=zeros(ndim,1); Input.v0=zeros(ndim,1); Input.a0=zeros(ndim,1);

%Damping parameters
%Rayleigh Damping with damping ratios (zetas) for the first two modes (OmegaIndexes)
Input.zeta = [0.02 0.02]; Input.OmegaIndexes = [1 2];
% See: GetRayleighDamping function to change this
%If damping is left in comment the damping parameters of the Input file are
%implemented!!!!!!
%To directly set Rayleigh parameters modify the input files (last two lines)
% and comment out the zeta,OmegaIndexes definition


% Time integration parameters - Notation similar to SDOF benchmark (see description)
fs = 100;              % working sampling frequency.
upsamp = 1;            % upsampling factor to ensure integration accuracy.
fsint = fs*upsamp;     % integration sampling frequency.
Input.dt = 1/fsint;    % integration time step.

%% Excitation signal design. 

%Example inputs are provided for demonstration along with a multisine constructor
method = 'ExampleB';   % Alternatives 'ExampleA' / 'ExampleB' / 'sinus'
Angle=pi/4;            % Angle of motion. See ReadMe documentation for alternative definitions
AmpF =0.1;            % Amplitude coefficient.
P = 3;                 % number of excitation periods.
N = 500;               % number of points per period.

fmin = 3;              % excitation bandwidth.
fmax = 50;

[SynthesizedAccelerogram,lowpass] = ExcitationDesign(method,...
    AmpF, upsamp, P, N, fmin, fmax,fsint);
Input.dt = 1/fsint;

%So, by changing this the user can speficy a different excitation (ground motion acceleration signal)
Input.SynthesizedAccelerogram=SynthesizedAccelerogram; 
Input.Angle=Angle;


%% Evaluate model
MODEL = BoucWenRun(Input, MODEL);

%% Low-pass filtering and downsampling in case upsampling was employed
% Software follows the notation in SDOF benchmark (see description)

if upsamp>1
    MODELV.Uups = MODELV.U;
    MODELV.U = Downsampling(MODELV.U,upsamp,P,N);
    MODEL.Uups = MODEL.U;
    MODEL.U = Downsampling(MODEL.U,upsamp,P,N);
    f = downsample(f,upsamp);
    % Removal of the last simulated period to eliminate the edge effects
    %due to the low-pass filter.
    f = f(1:(P-1)*N,:);
    Input.SynthesizedAccelerogramUps=f;
    Input.SynthesizedAccelerogram=f;
    P = P-1;
end

%% Results visualization
% In postprocessing you need to check the proper dof to visualize the 
% hysteresis curve as not all spring dofs are activated

%Plot deformed state of frame
%MODEL.u is the plotted deformed configuration. 
%Currently includes last timestep but can be freely modified.
% scale = 3*1/max(abs(MODEL.u([1:6:end,2:6:end,3:6:end])));
% plot_model( MODELC, scale )

%Plot time history
[~,ndof]=max(max(abs(MODEL.U),[],2));
plot(MODEL.U(ndof,:));

%Plot hysterisis loop on one of the nonlinear links
% [~,ndof]=max(max(abs(MODELC.HistR),[],2));
ndof=5;
figure
plot(MODEL.HistU(ndof, 1:(end-1)),MODEL.HistR(ndof,1:(end-1)),'-b');