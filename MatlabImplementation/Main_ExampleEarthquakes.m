%% Main file to evaluate the response of benchmark simulator
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

%% Clear workspace and load path 
clc;clear;
addpath InputFiles
addpath core

%% Define Input
% Definition of input files for the analysis
% MODEL = InputFileLinksBeams();    
% MODEL = InputFileLinksFirst_DifferentBW();    
% MODEL = InputFileLinksFirst();    
% MODEL = InputFileLinksAll();
% MODEL= AutomaticMesh(3,2,2,[5.0, 4.0 ,3.0]);

% Plot undeformed state of the model and visualize beam orientation
plot_model( MODEL, 0 )

% Initial conditions (displacements/velocities)
Input.u0=zeros(MODEL.ndim,1); Input.v0=zeros(MODEL.ndim,1); 

% Damping parameters (MODEL.dyn.a, MODEL.dyn.b)
% Default: Rayleigh damping with coefficients defined in Input file.
% Customized: Define damping ratios (zetas) & mode number (OmegaIndexes)
MODEL.zeta = [0.02 0.02]; MODEL.OmegaIndexes = [1 2]; 

%% Define Bouc-Wen parameters

%Bouc-Wen hysteresis model parameters - Notation follows description pdf
%bw_a=a|Alpha=A|N=w|Beta=b|Gamma=g|deltav=dv|deltan=dnu|bw_k=k

BW.bw_a = 0.25; BW.Alpha=1.0; BW.N=1; BW.Beta=3; BW.Gamma=2; 
BW.deltav = 0; BW.deltan=0; BW.bw_k = 2.0e08;
BW.integration_method='RK4'; % 'Euler' | 'RK2' | 'RK4'
MODEL.BW=BW;

%% Excitation

% Time integration parameters
fintegration = 1000;          % working sampling frequency.

%Excitation generator parameters
method = 'Morgan';             % Alternatives 'Kobe'|'ElCentro'|'Morgan'
AmpF =1e6;                     % Amplitude coefficient.
angle=pi/4;                  % Angle of ground motion.

SynthesizedAccelerogram = ExcitationDesignExamples(method,AmpF,fintegration);

%So, by changing this the user can speficy a different excitation (ground motion acceleration signal)
Input.SynthesizedAccelerogram=SynthesizedAccelerogram';
Input.angle=angle;
MODEL.dyn.nt = length(Input.SynthesizedAccelerogram); 
MODEL.dyn.dt = 1/fintegration;
Input.dt = MODEL.dyn.dt; t = 0:MODEL.dyn.dt:(MODEL.dyn.nt-1)*MODEL.dyn.dt;

%% Evaluate model
MODEL = Evaluate_Frame(MODEL, Input);

%% Results visualization

%Plot time history
[~,ndof]=max(max(abs(MODEL.U),[],2));
figure
plot(t, MODEL.U(ndof,:));

%Plot time history
[~,ndof]=max(max(abs(MODEL.Uorig),[],2));
figure
plot(t, MODEL.Uorig(ndof,:));

%Plot hysterisis loop on one of the nonlinear links
[~,ndof]=max(max(abs(MODEL.BW.HistR),[],2));
figure
plot(MODEL.BW.HistU(ndof, 1:(end-1)),MODEL.BW.HistR(ndof,1:(end-1)),'-b');

%Plot phase space
[~,ndof]=max(max(abs(MODEL.U),[],2));
figure
plot(MODEL.U(ndof,:),MODEL.V(ndof,:))

%Plot deformed state of frame
%MODEL.u is the plotted deformed configuration. 
%Currently includes last timestep but can be freely modified.
% scale = 3*1/max(abs(MODEL.u([1:6:end,2:6:end,3:6:end])));
% plot_model( MODEL, scale )