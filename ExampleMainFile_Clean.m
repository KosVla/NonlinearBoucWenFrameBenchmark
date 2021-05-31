%% Example main file to evaluate the response of two-story shear frame 
% with BW links and demonstrate all possible inputs
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.
%and
% J. Noel and M. Schoukens, 
% Hysteretic benchmark with a dynamic nonlinearity,
% in Workshop on non-linear system identification benchmarks, 2016, pp. 7–14

%% Clear workspace and load path 
clc;clear;
addpath InputFiles
addpath core

%% Define Input
% Definition of input files for the analysis

%Select Model 

% MODEL = InputFileLinksBeams();                %Model with only horizontal beam links
% MODEL = InputFileLinksFirst();                %Additional links on the basement columns
% MODEL = InputFileLinksFirst_DifferentBW();    %Model with different BC parameters for a link
MODEL = InputFileLinksAll();                    %Links on all beams and columns

ndim=6*numel(MODEL.nodes(:,1)); 
%Plot undeformed state of the model and visualize beam orientation
% plot_model( MODEL, 0 )


%Define parametric dependencies

%Bouc-Wen hysteresis model parameters - Notation follows description pdf
%bw_a=a / Alpha=A / N=w / Beta=b / Gamma=g / deltav=dv / deltan=dnu
%bw_k=k
Input.bw_a = 0.10; Input.Alpha=1.0;
Input.N=1; Input.Beta=3; Input.Gamma=2; Input.deltav = 0; Input.deltan=0;
Input.bw_k = 1.62e08;
%Additional amplitude parameter for BW. Multiplies only hysteretic term to
%magnify influence
Input.AmpBW=1.0;

%Initial conditions (displacements/velocities/accelerations)
Input.u0=zeros(ndim,1); Input.v0=zeros(ndim,1); Input.a0=zeros(ndim,1);

%Damping parameters
%Rayleigh Damping with damping ratios (zetas) for the first two modes (OmegaIndexes)
Input.zeta = [0.01 0.01]; Input.OmegaIndexes = [1 2]; 
% See: GetRayleighDamping function to change this
%If damping is left in comment the damping parameters of the Input file are
%implemented!!!!!!
%To directly set Rayleigh parameters modify the input files (last two lines)
% and comment out the zeta,OmegaIndexes definition

% Time integration parameters - Notation similar to SDOF benchmark (see description)
fs = 100;              % working sampling frequency.
upsamp = 2;            % upsampling factor to ensure integration accuracy.
fsint = fs*upsamp;     % integration sampling frequency.

%% Excitation signal design. 

%Example inputs are provided for demonstration along with a multisine constructor
method = 'sinus';   % Alternatives 'ExampleA' / 'ExampleB' / 'sinus'
Angle=pi/4;         % Angle of motion. See ReadMe documentation for alternative definitions
AmpF =1e5;            % Amplitude coefficient.
P = 1;              % number of excitation periods.
N = 1000;           % number of points per period.

fmin = 3;           % excitation bandwidth.
fmax = 50;

[SynthesizedAccelerogram,lowpass] = ExcitationDesign(method,...
    AmpF, upsamp, P, N, fmin, fmax,fsint);
Input.dt = 1/fsint;

%So, by changing this the user can speficy a different excitation (ground motion acceleration signal)
Input.SynthesizedAccelerogram=SynthesizedAccelerogram; 
Input.Angle=Angle;


%% Evaluate model
[MODEL] = BoucWenRun(Input,MODEL);

%% Low-pass filtering and downsampling in case upsampling was employed
% Software follows the notation in SDOF benchmark (see description)

if lowpass==1
    MODEL.Uups = MODEL.U;
    MODEL.U = Downsampling(MODEL.U,upsamp,P+1,N);
    MODEL.Vups = MODEL.V; MODEL.Aups = MODEL.A;
    MODEL.A = Downsampling(MODEL.A,upsamp,P+1,N);
    MODEL.V = Downsampling(MODEL.V,upsamp,P+1,N);
    f = downsample(Input.SynthesizedAccelerogram,upsamp);
    % Removal of the last simulated period to eliminate the edge effects
    %due to the low-pass filter.
    f = f(1:(P+1-1)*N,:);
    Input.SynthesizedAccelerogramUps=f;
    Input.SynthesizedAccelerogram=f;
    P = P+1-1;
end

%% Results visualization

%Plot deformed state of frame
%MODEL.u is the plotted deformed configuration. 
%Currently includes last timestep but can be freely modified.
% scale = 3*1/max(abs(MODEL.u([1:6:end,2:6:end,3:6:end])));
% plot_model( MODEL, scale )

%Plot time history
ndof=1;
plot(MODEL.U(ndof,:));


%Plot hysterisis loop on one of the nonlinear links
ndof=1;
figure
plot(MODEL.HistU(ndof, 1:(end-1)),MODEL.HistR(ndof,1:(end-1)),'-b');
