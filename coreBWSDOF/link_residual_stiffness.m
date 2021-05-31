function [Re, Ke, link_hist] = link_residual_stiffness(ue, link_properties, link_hist)
%Evaluates the nonlinear link element stiffness matrix and internal force vector
%
%Inputs:
%  ue : vector / Displacements of the current timestep on the dofs of 
%               the element
%  link_properties : vector / Properties of the Bouc-Wen model of the link
%  link_hist : struct / The Bouc-Wen history parameters-terms of the previous step
%
%Returns:
%  Re : vector / Updated restoring forces of the hysteretic link
%  Ke : matrix / Stiffness matrix of the hysteretic link
%  link_hist: struct / Updated history struct for the Bouc-Wen model
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.


Ke = zeros(12);
Re = zeros(12,1);

Ki = [1 -1; -1 1];
Ri = [1; -1];
       
for d=1:6
    du = ue(d)-ue(d+6);
    d_prop = link_properties(d,:);

    prop.Amp = d_prop(7); prop.a = d_prop(1);
    prop.k = d_prop(2); prop.Alpha=d_prop(3);
    prop.Beta=d_prop(4); prop.Gamma=d_prop(5);
    prop.N=d_prop(6); prop.deltav=d_prop(8);
    prop.deltan=d_prop(9);
    
    [ rbw,kbw,link_hist(d) ] = bw_1d( du, link_hist(d), prop);
        
    Ke([d d+6],[d d+6])=Ke([d d+6],[d d+6])+kbw*Ki;
    Re([d d+6])=Re([d d+6]) + rbw*Ri;
    
end

end