function [alpha,beta,C] = GetRayleighDamping(Kstiff,Mmass,zeta,OmegaIndexes)
%Formulate Rayleigh Damping for the system based on damping ratios 
% and modal frequencies contribution
%
%Inputs:
% Kstiff,Mmass: matrices / Stiffness and mass matrices of the system
% zeta : vector / Vector containing the damping ratios for each frequency
% OmegaIndexes : vector / Vector containing the indexes of the Omegas
%                           to apply the damping ratios
%Returns:
% alpha,beta: floats / Rayleigh damping coefficients. 
%           Alpha for mass proportionality, Beta for stiffness proportionality.
% C: matrix / Damping matrix based on the Rayleigh formulation
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

if nargin==3
    no_modes=10;
else
    no_modes=max(OmegaIndexes)+2;
end

[~,values] = eigs(Kstiff,Mmass,no_modes,'SM');
Eig_values_mat=sort(diag(values));
%Get natural frequencies
Omegas = sqrt(Eig_values_mat);
% Natural_freq_mat = Omegas/2/pi;

Omega1 = Omegas(OmegaIndexes(1));
Omega2 = Omegas(OmegaIndexes(2));

if length(zeta)==1
    beta = 0;
    alpha = 2*Omegas(OmegaIndexes(1))*zeta(1);
elseif length(zeta)==2
    beta= (2*Omega2*zeta(2)-2*Omega1*zeta(1))/(Omega2^2-Omega1^2);
    alpha=2*zeta(1)*Omega1-beta*Omega1^2;
end

C = alpha *  Mmass + beta * Kstiff;

end