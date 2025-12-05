function [alpha,beta] = GetRayleighDamping(Kstiff,Mmass,zeta,OmegaIndexes)
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

no_modes=1;
condition = true;
while condition
    no_modes=no_modes*10;
    [~,values] = eigs(Kstiff,Mmass,no_modes,'SM');
    Eig_values_mat=sort(diag(values));
    %Get natural frequencies
    Omegas = sqrt(Eig_values_mat);
    condition = imag(Omegas(end))>0 || real(Omegas(end))<1e-3;
end
    
count=0;
condition = true;
while condition
    count=count+1;
    condition = imag(Omegas(count))>0 || real(Omegas(count))<1e-3;
end

if count==1
    Omega1 = Omegas(OmegaIndexes(1));
    Omega2 = Omegas(OmegaIndexes(2));
else
    Omega1 = Omegas(count);
    Omega2 = Omegas(count+1);
end

if isscalar(zeta)
    beta = 0;
    alpha = 2*Omegas(OmegaIndexes(1))*zeta(1);
elseif length(zeta)==2
    beta= (2*Omega2*zeta(2)-2*Omega1*zeta(1))/(Omega2^2-Omega1^2);
    alpha=2*zeta(1)*Omega1-beta*Omega1^2;
end

% C = alpha *  Mmass + beta * Kstiff;

end