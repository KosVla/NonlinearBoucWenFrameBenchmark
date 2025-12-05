function [ dofs ] = element_global_dofs( nodes )
%Assembles the global degrees of freedom of an element from the nodes
%
%Input:
%  nodes : vector / Nodes of the element (start and end node numbers)
%Returns:
%  dofs : vector / Degrees of freedom in the global reference system
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.


dofs = [];

for n=1:numel(nodes)
    dofs = [dofs (1+6*(nodes(n)-1)):(6+6*(nodes(n)-1))];
end

end

