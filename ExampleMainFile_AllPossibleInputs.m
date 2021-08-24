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
addpath core
addpath InputFiles

%% Define Input
% Definition of input files for the analysis

%Select Model 
% MODEL = InputFileLinksBeams();                %Model with only horizontal beam links
% MODEL = InputFileLinksFirst();                %Additional links on the basement columns
% MODEL = InputFileLinksFirst_DifferentBW();    %Model with different BC parameters for a link
% MODEL = InputFileLinksAll();                    %Links on all beams and columns

%Automatic Mesh examples
%NoofFloors excludes the basement, so floors additional to the basement

% Example structure same as InputFileLinksAll, node numbering is different!
NoOfFloors=1; NoOfFramesx=2; NoOfFramesy=1; dimensions=[7.50,5.0,3.2];
MODEL = AutomaticMesh(NoOfFloors,NoOfFramesx,NoOfFramesy,dimensions); 

% Example structure 2
% MODEL = AutomaticMesh(4,2,2,[7.50,5.0,3.2]); 

%Plot undeformed state of the model and visualize beam orientation
% plot_model( MODEL, 0 )

ndim=6*numel(MODEL.nodes(:,1)); 

%Define parametric dependencies

%Bouc-Wen hysteresis model parameters - Notation follows description pdf
%bw_a=a / Alpha=A / N=w / Beta=b / Gamma=g / deltav=dv / deltan=dnu
%bw_k=k
Input.bw_a = 0.10; Input.Alpha=1.0;
Input.N=1; Input.Beta=3; Input.Gamma=2; Input.deltav = 0; Input.deltan=0;
Input.bw_k = 1.62e6;
%Additional amplitude parameter for BW. Multiplies only hysteretic term to
%magnify influence
Input.AmpBW=1.0;

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

%Define excitation signal/ground motion acceleration and angle.
%Two example inputs are provided for demonstration along with a multisine
%constructor
method = 'ExampleB'; % Alternatives 'ExampleA' / 'ExampleB' / 'sinus'

%If the user wants to define the input signal manually edit the following:
%Input.SynthesizedAccelerogram = Input signal time history
%Input.Angle = Angle of excitation, see below for alternative definitions

%Nodal excitation can be defined by modifying the input file!!!

switch method
    case 'ExampleA'
        load('ExampleA')
        
        AmpF =10;               % Amplitude coefficient.
        Input.SynthesizedAccelerogram = AmpF*ExciteA;
        
        %Angle definition - Alternative possibilities
        % 2D angle in xy plane with axis x
        Input.Angle=pi/4;
        % 3D angle with the following notation:
        %[Angle(3D vector with z axis) Angle(XY plane projection with x Quadrant]
        % Input.Angle=[ pi/2 pi/4 1]; %Equivalent to the first alternative.
        % 3D vector [Coordx Coordy Coordz 0]
%         Input.Angle=[cos(pi/4) sin(pi/4) 0 0]; %Equivalent to the first
%         alternative.

    case 'ExampleB'
        load('ExampleB')
        AmpF =0.1;
        Input.SynthesizedAccelerogram = AmpF*ExciteB;
        Input.Angle=pi/4; 
    otherwise
        %Multisine constructor - Notation similar to SDOF benchmark (see description)

        lowpass=1;
        P = 3;                  % number of excitation periods.
        N = 1000;               % number of points per period.
        Nint = N*upsamp;        % number of points per period during integration.

        fmin = 3;               % excitation bandwidth.
        fmax = 10;
        
        A = 1e5;                 % excitation amplitude coefficient.
        

        Pfilter = 1;            % extra period to avoid edge effects during low-pass filtering.
        P = P + Pfilter;

        F = zeros(Nint,1);      % definition of the multisine excitation.
        fres = fsint/Nint;
        exclines = 1:ceil(fmax/fres);
        exclines(exclines < floor(fmin/fres)) = [];

        F(exclines+1) = exp(complex(0,2*pi*rand(size(exclines))));
        f = 2*real(ifft(F));
        f = A*f/std(f);
        f = repmat(f,[P 1]);
        Input.SynthesizedAccelerogram=f; 
        %Same alternatives as before apply for the angle
        Input.Angle=pi/4;
%         Input.Angle=[ pi/2 pi/4 1];
end

%% Evaluate model
[MODEL] = BoucWenRun(Input,MODEL);
Input.SynthesizedAccelerogram(end+1)=0;
%% Low-pass filtering and downsampling in case upsampling was employed
% Software follows the notation in SDOF benchmark (see description)

if upsamp>1
    MODEL.Uups = MODEL.U;
    MODEL.U = Downsampling(MODEL.U,upsamp,P,N);
    MODEL.Vups = MODEL.V; MODEL.Aups = MODEL.A;
    MODEL.A = Downsampling(MODEL.A,upsamp,P,N);
    MODEL.V = Downsampling(MODEL.V,upsamp,P,N);
    f = downsample(f,upsamp);
    % Removal of the last simulated period to eliminate the edge effects
    %due to the low-pass filter.
    f = f(1:(P-1)*N,:);
    Input.SynthesizedAccelerogramUps=f;
    Input.SynthesizedAccelerogram=f;
    P = P-1;
end

%% Results visualization

%Plot deformed state of frame
%MODEL.u is the plotted deformed configuration. 
%Currently includes last timestep but can be freely modified.
% scale = 3*1/max(abs(MODEL.u([1:6:end,2:6:end,3:6:end])));
% plot_model( MODEL, scale )

%Plot time history
[~,ndof]=max(max(abs(MODEL.U),[],2));
plot(MODEL.U(ndof,:));


%Plot hysterisis loop on one of the nonlinear links
[~,ndof]=max(max(abs(MODEL.HistR),[],2));
figure
plot(MODEL.HistU(ndof, 1:(end-1)),MODEL.HistR(ndof,1:(end-1)),'-b');
