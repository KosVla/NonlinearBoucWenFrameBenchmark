function [ MODEL ] = assemble_nlBW( MODEL )
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

nt = MODEL.nt;
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

% Assembly of the plate elements

if (isfield(MODEL,'plate_elements'))
    n_plates = numel(MODEL.plate_elements(:,1));
else
    n_plates = 0;
end

for e=1:n_plates
    plate_nodes = MODEL.plate_elements(e,1:end);
    plate_COORDS = MODEL.nodes(plate_nodes,:);
    plate_material_properties = MODEL.material_properties(MODEL.plate_material_properties(e),:);
    plate_thickness = MODEL.plate_thickness(e,:);
%     plate_loads = MODEL.plate_loads(e,:);

    %Check to avoid singularities
    if plate_COORDS(1,3)==plate_COORDS(2,3) && plate_COORDS(1,3)==plate_COORDS(4,3)
        plate_COORDS(3,3)= plate_COORDS(3,3)*0.9995;
        plate_COORDS(4,3)= plate_COORDS(4,3)*0.9995;
    end
    
    [Ke,Me] = shell_mass_stiffness_rhs( plate_COORDS, plate_material_properties, plate_thickness);
    [ plate_dofs ] = element_global_dofs( plate_nodes );
    ue = MODEL.u(plate_dofs);

    MODEL.K(plate_dofs,plate_dofs) = MODEL.K(plate_dofs,plate_dofs) + Ke;
    MODEL.M(plate_dofs,plate_dofs) = MODEL.M(plate_dofs,plate_dofs) + Me;
    MODEL.fint(plate_dofs) = MODEL.fint(plate_dofs) + Ke*ue;
%     MODEL.f(plate_dofs) = MODEL.f(plate_dofs) + fe;
end

% Assembly of nonlinear links
Properties=MODEL.DofsLoop; HIST=MODEL.HistBW;

DofsLoop=Properties(:,1:2);
% link_properties = MODEL.nl_link_bw_properties(1,:);
link_properties = Properties(:,3:end);

[ Rs,Ks, HIST ] = bw_Ndall( MODEL.u, DofsLoop, HIST, link_properties );
MODEL.HistBWtemp = HIST;
MODEL.HistR(:,nt)=Rs; MODEL.HistU(:,nt)=HIST.Um;

Fnl =sparse( [DofsLoop(:,1); DofsLoop(:,2)], ones(size(DofsLoop,1)*2,1), [Rs; -Rs], ndof, 1);
MODEL.fint=MODEL.fint+Fnl;

rows = [ DofsLoop(:,1); DofsLoop(:,1); DofsLoop(:,2); DofsLoop(:,2) ];
columns = [ DofsLoop(:,1); DofsLoop(:,2); DofsLoop(:,1); DofsLoop(:,2) ];
values=[ Ks; -Ks; -Ks; Ks];
Knl =sparse( rows, columns, values, ndof, ndof);

MODEL.K = MODEL.K + Knl;

%Assembly of degrees of freedom not activated
% if (isfield(MODEL,'DofsLoopLinear'))
%     DofsLoopL=MODEL.DofsLoopLinear;
%     
%     %Employ linear link
%     link_properties(1)=1.0; link_properties(7)=0.0;
%     [ rls,ks, MODEL.HistBWL ] = bw_Ndall( MODEL.u, DofsLoopL, MODEL.HistBWL, link_properties );
%     
%     rows = [ DofsLoopL(:,1); DofsLoopL(:,1); DofsLoopL(:,2); DofsLoopL(:,2) ];
%     columns = [ DofsLoopL(:,1); DofsLoopL(:,2); DofsLoopL(:,1); DofsLoopL(:,2) ];
%     values=[ ks; -ks; -ks; ks];
%     Kl =sparse( rows, columns, values, ndof, ndof);
%     MODEL.K = MODEL.K + Kl;
%     values=[ rls; -rls];
%     Fl =sparse( [DofsLoopL(:,1); DofsLoopL(:,2)], ones(size(DofsLoopL,1)*2,1), values, ndof, 1);
%     MODEL.fint=MODEL.fint+Fl;
% end


% Assembly of the springs

if (isfield(MODEL,'springs')&&(numel(MODEL.springs)>0))

    n_springs = numel(MODEL.springs(:,1));
else
    n_springs = 0;
end

for s=1:n_springs
    node = MODEL.springs(s,1);
    spring_constants = MODEL.springs(s,2:end);
    [ dofs ] = node_global_dofs( node );
    MODEL.K(dofs,dofs) = MODEL.K(dofs,dofs) + diag(eye(numel(dofs))*spring_constants');
end

% Assembly of masses

if (isfield(MODEL,'masses')&&(numel(MODEL.masses)>0))

    n_masses = numel(MODEL.masses(:,1));
else
    n_masses = 0;
end

for m=1:n_masses
    node = MODEL.masses(m,1);
    mass = MODEL.masses(m,2:end);
    [ dofs ] = node_global_dofs( node );
    MODEL.M(dofs,dofs) = MODEL.M(dofs,dofs) + diag(eye(numel(dofs))*mass');
end

end
