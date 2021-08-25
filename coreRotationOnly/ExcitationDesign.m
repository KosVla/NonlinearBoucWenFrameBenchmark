function [SynthesizedAccelerogram,lowpass] = ExcitationDesign(...
    method,AmpF, P, N,upsamp, fmin, fmax,fsint)
%Function to assemble the excitation ground acceleration signal 
%
%Parameters:
%   method : string / Define type of signal/ground motion acceleration.
%                     Two example inputs are provided for demonstration
%                     along with a multisine constructor. 
%                     Alternatives 'ExampleA' / 'ExampleB' / 'sinus'
%   AmpF : floats / Amplitude coefficient for signal 
%   P,N,upsamp,fmin,fmax, fsint : floats / Parameters for the multisine.
%                            Number of excitaion periods, number of points 
%                            per period, upsampling factor, minimum, 
%                            maximum and sampling frequency of signal
%
%
%Returns:
%   SynthesizedAccelerogram : vector / Contains the synthesized signal 
%
%
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.
%and
% J. Noel and M. Schoukens, 
% Hysteretic benchmark with a dynamic nonlinearity,
% in Workshop on non-linear system identification benchmarks, 2016, pp. 7–14

%Variable to control if lowpass filtering is needed after the integration
%due to upsampling 
lowpass=0;


%If the user wants to define the input signal manually edit the following:
%Input.SynthesizedAccelerogram = Input signal time history
%Input.Angle = Angle of excitation, see below for alternative definitions

%Nodal excitation can be defined by modifying the input file!!!

switch method
    case 'ExampleA'
        load('ExampleA')
        SynthesizedAccelerogram = AmpF*ExciteA;
    case 'ExampleB'
        load('ExampleB')
        SynthesizedAccelerogram = AmpF*ExciteB; 
    otherwise
        %Multisine constructor - Notation similar to SDOF benchmark (see description)

        lowpass=1;
        Nint = N*upsamp;        % number of points per period during integration.
        
        A = AmpF;                 % excitation amplitude coefficient.
        

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
        SynthesizedAccelerogram=f; 
end

end

