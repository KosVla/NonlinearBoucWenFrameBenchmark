function [Error] = CheckErrorStruct(MODEL,MODELR)
%Function to check error between models
%
%Returns:
%  Error : struct / The Error struct contains the error in displacements 
%                   (Error.Norm), in velocities (Error.NormV) and in
%                   (nonlinear) restoring forces (Error.NormR)
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

ndof = size(MODEL.U,1);

U = MODEL.U(1:ndof,1:end-1); Ur=MODELR.U(1:ndof,1:end-1);
V = MODEL.V(1:ndof,1:end-1); Vr=MODELR.V(1:ndof,1:end-1);

Error.Norm = norm(U(:)-Ur(:))/norm(U(:));
Error.NormV = norm(V(:)-Vr(:))/norm(V(:));

NLforces = MODEL.HistR; 
NLforcesR = MODELR.HistR; 

Error.NormR = norm(NLforces(:)-NLforcesR(:))/norm(NLforces(:));
end

