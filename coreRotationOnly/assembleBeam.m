function [ MODEL ] = assembleBeam( MODEL )
%Assembles the updated system matrices based on the displacement of the current timestep
%
%Input:
%  MODEL : struct / The MODEL struct contains all system properties, 
%                   parameters and matrices. Here, the MODEL.u represents
%                   the displacement of the current timestep and is used
%                   as an input variable to update the internal forces and
%                   compute the Bouc-Wen hysteretic forces of the links. 
%                   The Bouc-Wen variables of the previous step are located
%                   in MODEL.HistBW. Additional properties stored in the 
%                   MODEL struct are also used.
%
%Returns:
%  MODEL : struct / The MODEL struct contains all system properties,
%                   parameters and matrices. The stiffness matrix and the 
%                   internal forces vector are updated, along with the
%       			respective history variables of the Bouc-Wen model
%                   (zeta, x etc -> MODEL.HistBW)
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

ndof = 6*numel(MODEL.nodes(:,1));

MODEL.K = sparse(ndof,ndof);
MODEL.M = sparse(ndof,ndof);
MODEL.fint = zeros(ndof,1);
MODEL.f = zeros(ndof,1);

condition=false;
if (isfield(MODEL,'u')==0)
    MODEL.u = zeros(ndof,1);
    condition=true;
end

% Assembly of the beam elements

n_beams =size(MODEL.beam_elements,1);

for e=1:n_beams
    if condition
        %beam_type = MODEL.beam_elements(e,1);
        beam_nodes = MODEL.beam_elements(e,2:end);
        beam_COORDS = MODEL.nodes(beam_nodes,:);
        beam_material_properties = MODEL.material_properties(MODEL.beam_material_properties(e),:);
        beam_cross_section_properties = MODEL.cross_sections(MODEL.beam_cross_sections(e,:)',:);
        beam_loads = MODEL.beam_loads(e,:);

        [ beam_dofs ] = element_global_dofs( beam_nodes );
        [Ke,Me,fe] = beam_mass_stiffness_rhs(beam_COORDS, beam_material_properties, beam_cross_section_properties, beam_loads);
        MODEL.ElementsK{e}=Ke;
        MODEL.ElementsDofs{e}=beam_dofs;
        MODEL.K(beam_dofs,beam_dofs) = MODEL.K(beam_dofs,beam_dofs) + Ke;
        MODEL.M(beam_dofs,beam_dofs) = MODEL.M(beam_dofs,beam_dofs) + Me;
        MODEL.f(beam_dofs) = MODEL.f(beam_dofs) + fe;
    else
        beam_dofs=MODEL.ElementsDofs{e};
        Ke = MODEL.ElementsK{e};
        MODEL.K(beam_dofs,beam_dofs) = MODEL.K(beam_dofs,beam_dofs) + Ke;
        MODEL.M=MODEL.Mbeam;
        MODEL.f=MODEL.fbeam;
    end
    
    ue = MODEL.u(beam_dofs);
    MODEL.fint(beam_dofs) = MODEL.fint(beam_dofs) + Ke*ue;
end

if condition
    MODEL.Kbeam=MODEL.K;
    MODEL.Mbeam=MODEL.M;
    MODEL.fbeam=MODEL.f;
end


end
