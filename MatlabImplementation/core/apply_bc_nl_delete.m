function [ MODEL ] = apply_bc_nl_delete( MODEL )
%Applies the boundary conditions onthe system matrices
%
%Input:
%  MODEL : struct / The MODEL struct contains all system properties, parameters and matrices
%			The function applies the boundary conditions on the stiffness, mass and forces matrices
%
%Returns:
%  MODEL : struct / The MODEL struct contains all system properties, parameters and matrices
%			The respective matrices are updated and stored in the same entries. 
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

ndofs = numel(MODEL.nodes(:,1))*6;

% Nodal loads
if (numel(MODEL.nodal_loads)>0)

    n_loads = numel(MODEL.nodal_loads(:,1));

    for l=1:n_loads
        node = MODEL.nodal_loads(l,1);
        load = MODEL.nodal_loads(l,2:end);
        [ dofs ] = node_global_dofs( node );
        MODEL.f(dofs) = MODEL.f(dofs) + load';
    end

end

% Imposed displacements

if (numel(MODEL.nodal_displacements)>0)
    
    bc_dofs = [];
    bc_dofs_nz = [];
    bc_values = [];

    for d=1:numel(MODEL.nodal_displacements(:,1))
        node = MODEL.nodal_displacements(d,1);
        dofs= node_global_dofs( node );
        dofs_act = dofs(MODEL.nodal_displacements(d,2:7)~=0);        
        nz = MODEL.nodal_displacements(d,8:end)~=0 ;
        nz(nz==0)=[];
        if (numel(nz)>0)
            bc_dofs_nz = [bc_dofs_nz; dofs(nz)'];
            bc_values = [bc_values; MODEL.nodal_displacements(d,nz+7)'];
            dofs_act(dofs_act==dofs(nz))=[]; 
        end       
        bc_dofs  = [bc_dofs; dofs_act'];        
    end
    
    if (numel(bc_dofs_nz)>0)
        fbc = MODEL.K(:,bc_dofs_nz)*bc_values;
        MODEL.f = MODEL.f - fbc;
        MODEL.fint(bc_dofs_nz) = 0;
        MODEL.u(bc_dofs_nz)=bc_values;
    end
    
    MODEL.K(bc_dofs, :) = [];
    MODEL.K(:, bc_dofs) = [];
    MODEL.M(bc_dofs, :) = [];
    MODEL.M(:, bc_dofs) = [];
    MODEL.f(bc_dofs) = [];
    MODEL.fint(bc_dofs) = [];
%     MODEL.u(bc_dofs) = [];
    MODEL.freedofsN = 1:ndofs; 
    MODEL.freedofsN(bc_dofs)=[];
    MODEL.freedofs=1:size(MODEL.K,1);
end

end