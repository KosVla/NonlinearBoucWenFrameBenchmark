function [ dofs ] = node_global_dofs( node )
%Evaluates the degrees of freedom from the node number
%
%Inputs:
% node: float / Node number on the respective coordinates matrix
%
%Returns:
% dofs: vector / Numbering of degrees of freedom on the input node
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.


dofs = (1+6*(node-1)):(6+6*(node-1));

end