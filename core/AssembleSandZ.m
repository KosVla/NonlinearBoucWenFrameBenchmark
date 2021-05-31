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
if length(MODEL.angle)==1 
    % Assumes input notation [Angle in xy plane with x axis]
    Zacceler=InputParams.SynthesizedAccelerogram;
    originalnodes = size(MODEL.nodes,1)-size(MODEL.nl_link_elements,1);
    loaddofsx = ([1:originalnodes]-1)*6+1;
    loaddofsy=([1:originalnodes]-1)*6+2;
    MODEL.Loaddofs = reshape([loaddofsx; loaddofsy],originalnodes*2,1);
    Mg = diag(GetLumpedMass(MODEL.M, [1:6:neq 2:6:neq 3:6:neq]));
    % Mg=diag(MODEL.M); 
    S=sparse(length(Mg),1); 
    S(MODEL.Loaddofs)= Mg(MODEL.Loaddofs);
    S(MODEL.Loaddofs(1:2:end)) = S(MODEL.Loaddofs(1:2:end))*cos(MODEL.angle);
    S(MODEL.Loaddofs(2:2:end)) = S(MODEL.Loaddofs(2:2:end))*sin(MODEL.angle);
    F=sparse(Zacceler);
elseif length(MODEL.angle)==3
    % Assumes input notation
    % [Angle(3D vector with z axis) Angle(XY plane projection with x Quadrant]
    Zacceler=InputParams.SynthesizedAccelerogram;
    originalnodes = size(MODEL.nodes,1)-size(MODEL.nl_link_elements,1);
    loaddofsx = ([1:originalnodes]-1)*6+1;
    loaddofsy=([1:originalnodes]-1)*6+2;
    loaddofsz=([1:originalnodes]-1)*6+3;
    MODEL.Loaddofs = reshape([loaddofsx; loaddofsy; loaddofsz],originalnodes*3,1);
    Mg = diag(GetLumpedMass(MODEL.M, [1:6:neq 2:6:neq 3:6:neq]));
    % Mg=diag(MODEL.M); 
    S=sparse(length(Mg),1); 
    phi = MODEL.angle(1); theta=MODEL.angle(2);
    S(MODEL.Loaddofs)= Mg(MODEL.Loaddofs);
    switch MODEL.angle(3)
        case 1
            S(MODEL.Loaddofs(1:3:end)) = S(MODEL.Loaddofs(1:3:end))*sin(phi)*cos(theta);
            S(MODEL.Loaddofs(2:3:end)) = S(MODEL.Loaddofs(2:3:end))*sin(phi)*sin(theta);
            S(MODEL.Loaddofs(3:3:end)) = S(MODEL.Loaddofs(3:3:end))*cos(phi);
        case 2
            S(MODEL.Loaddofs(1:3:end)) = -S(MODEL.Loaddofs(1:3:end))*sin(phi)*cos(theta);
            S(MODEL.Loaddofs(2:3:end)) = S(MODEL.Loaddofs(2:3:end))*sin(phi)*sin(theta);
            S(MODEL.Loaddofs(3:3:end)) = S(MODEL.Loaddofs(3:3:end))*cos(phi);
        case 3
            S(MODEL.Loaddofs(1:3:end)) = -S(MODEL.Loaddofs(1:3:end))*sin(phi)*cos(theta);
            S(MODEL.Loaddofs(2:3:end)) = -S(MODEL.Loaddofs(2:3:end))*sin(phi)*sin(theta);
            S(MODEL.Loaddofs(3:3:end)) = S(MODEL.Loaddofs(3:3:end))*cos(phi);
        case 4
            S(MODEL.Loaddofs(1:3:end)) = S(MODEL.Loaddofs(1:3:end))*sin(phi)*cos(theta);
            S(MODEL.Loaddofs(2:3:end)) = -S(MODEL.Loaddofs(2:3:end))*sin(phi)*sin(theta);
            S(MODEL.Loaddofs(3:3:end)) = S(MODEL.Loaddofs(3:3:end))*cos(phi);
    end
    F=sparse(Zacceler);
elseif length(MODEL.angle)==4
    % Assumes input notation [Coordx Coordy Coordz 0] of vector of motion
    Coordx = MODEL.angle(1); Coordy=MODEL.angle(2); Coordz=MODEL.angle(3);
    Zacceler=InputParams.SynthesizedAccelerogram;
    originalnodes = size(MODEL.nodes,1)-size(MODEL.nl_link_elements,1);
    loaddofsx = ([1:originalnodes]-1)*6+1;
    loaddofsy=([1:originalnodes]-1)*6+2;
    loaddofsz=([1:originalnodes]-1)*6+3;
    MODEL.Loaddofs = reshape([loaddofsx; loaddofsy; loaddofsz],originalnodes*3,1);
    Mg = diag(GetLumpedMass(MODEL.M, [1:6:neq 2:6:neq 3:6:neq]));
    % Mg=diag(MODEL.M); 
    S=sparse(length(Mg),1); 
    S(MODEL.Loaddofs)= Mg(MODEL.Loaddofs);
    S(MODEL.Loaddofs(1:3:end)) = S(MODEL.Loaddofs(1:3:end))*Coordx;
    S(MODEL.Loaddofs(2:3:end)) = S(MODEL.Loaddofs(2:3:end))*Coordy;
    S(MODEL.Loaddofs(3:3:end)) = S(MODEL.Loaddofs(3:3:end))*Coordz;
    F=sparse(Zacceler);
end

end

