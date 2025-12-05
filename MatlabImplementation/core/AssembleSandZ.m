function [S,F] = AssembleSandZ(InputParams,MODEL)
%Assembles the excitation time history for all degrees of freedom - RHS of
%equation of motion
%
%Input:
%  MODEL : struct / The MODEL struct contains all system properties, 
%           parameters and matrices. Here the mass matrix is used to
%           assemble the lumped mass and then compute the nodal
%           contributions of the forcing. 
%  InputParams: struct / Input struct containing excitation relevant
%           parameters. Here the ground motion acceleration is used to
%           assemble the respective loading time history
%Returns:
%  F : matrix / Excitation time history
%  S : matrix / Weighting factor for each degree of freedom based on lumped
%               masses and angle of excitation. The RHS will be S*F'.
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

neq = numel(MODEL.u);
% Assumes input notation [Angle in xy plane with x axis]
Zacceler=InputParams.SynthesizedAccelerogram;
originalnodes = 1:(size(MODEL.nodes,1)- size(MODEL.nl_link_elements,1));

loaddofsx = (originalnodes-1)*6+1;
loaddofsy=(originalnodes-1)*6+2;
MODEL.Loaddofs = reshape([loaddofsx; loaddofsy],length(originalnodes)*2,1);
Mg = diag(GetLumpedMass(MODEL.Mall, [1:6:neq 2:6:neq 3:6:neq]));
S=sparse(length(Mg),1); 
S(MODEL.Loaddofs)= Mg(MODEL.Loaddofs);
S(MODEL.Loaddofs(1:2:end)) = S(MODEL.Loaddofs(1:2:end))*cos(InputParams.angle);
S(MODEL.Loaddofs(2:2:end)) = S(MODEL.Loaddofs(2:2:end))*sin(InputParams.angle);
F=sparse(Zacceler);

end

