function [SynthesizedAccelerogram] = ExcitationDesignExamples(...
    method,AmpF,fintegration)
% Assemble the ground acceleration signal used as excitation
%
%Parameters:
%   method : string / Define type of signal/ground motion acceleration.
%   AmpF : float / Amplitude coefficient for signal 
%   fintegration: float / Integration (sampling) frequency
%Returns:
%   SynthesizedAccelerogram : vector / Contains the synthesized signal 

switch method
    case 'Kobe'
        % load('KobeEarthq.mat')
        load("KobeAccelNoScaling.mat")
        Excitation = AmpF*KobeAccel(2,:); 
        fsample = KobeAccel(1,:);
    case 'ElCentro'
        load('ElCentroAccelNoScaling.mat')
        Excitation = AmpF*ElCentroAccel(2,:); 
        fsample = ElCentroAccel(1,:);
    case 'Morgan'
        load("MorganAccelNoScaling.mat")
        Excitation = AmpF*MorganAccel(2,:);
        fsample = MorganAccel(1,:);
end

xq = 0:1/fintegration:fsample(end);

if (fsample(2)<1/fintegration)
    fprintf('Warning: Integration timestep is larger than timestep of recordings')
end

SynthesizedAccelerogram = interp1(fsample,Excitation,xq);

end

