function [Mlumped] = GetLumpedMass(Mconsist, translational_dofs)
%Compute lumped mass matrix from consistent mass matrix based 
%on HRZ lumping scheme
%
%Literature:
%	Carlos Felipa, Introduction to Finite Element Methods, Chapter 31
%
%Inputs:
%   Mconsist: matrix / The consistent mass matrix
%   translational_dofs: matrix / The numbering of the dofs representing
%                                   translational degrees of freedom
%
%Returns:
%   Mlumped: matrix / The lumped mass matrix of the system
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.


%Extract diagonal coefficients of M consistent
diagonal_coeff = diag(Mconsist);

%Compute total mass
total_mass = sum(Mconsist,'all');

%Get number s of diagonal coefficients of translational dofs
coeff_s = sum(diagonal_coeff(translational_dofs));

%Compute diagonal of lumped matrix
lumped_diag = diagonal_coeff * total_mass/coeff_s;

%Compute Lumped Mass Matrix
Mlumped = diag(lumped_diag);

end

