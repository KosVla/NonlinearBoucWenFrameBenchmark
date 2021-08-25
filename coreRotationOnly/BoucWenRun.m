function [MODEL] = BoucWenRun(Input,MODEL)
%Main function assembling the input properly and integrating the MODEL
% forward in time
%
%Input:
%  Input : struct / Contains all input parameters relevant to Bouc-Wen, 
%                   time integration, Boundary conditions etc.
%  MODEL : struct / The MODEL struct contains all system properties,
%                   parameters and matrices. Here, the vectors and matrices
%                   needed for the integration are assembled and the 
%                   struct is given as input to the Newmark integration scheme
%Returns:
%  MODEL : struct / The MODEL struct contains all system properties, 
%                   parameters and matrices. Sfter integrating your model
%                   forward in time the MODEL struct contains the time
%                   histories of the displacement, velocities and 
%                   accelerations for every degree of freedom
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.


MODEL.angle= Input.Angle; MODEL.delete=true;
MODEL.dyn.dt = Input.dt; %time step
MODEL.dyn.nt = length(Input.SynthesizedAccelerogram); %duration of the analysis

MODEL = GetHystereticParamsRot(MODEL,Input);

tic;
[ MODEL ] = newmark_nlBW( MODEL );
time=toc;
MODEL.time=time;

if (isfield(MODEL,'nl_link_elements'))
    originaldofs = 6*(size(MODEL.nodes,1)-size(MODEL.nl_link_elements,1));
    MODEL.Uorig= MODEL.U(1:originaldofs,:);
else
    MODEL.Uorig= MODEL.U;
end

%Reconstruct results to include bounded degrees of freedom
if length(MODEL.u)~=6*numel(MODEL.nodes(:,1))
    temp=MODEL.u;
    MODEL.u = zeros(6*numel(MODEL.nodes(:,1)),1);
    MODEL.u(MODEL.freedofs) = temp;
    temp2=MODEL.U;
    MODEL.U = zeros(6*numel(MODEL.nodes(:,1)),size(temp2,2));
    MODEL.U(MODEL.freedofs,:) = temp2;
    temp3=MODEL.V;
    MODEL.V = zeros(6*numel(MODEL.nodes(:,1)),size(temp2,2));
    MODEL.V(MODEL.freedofs,:) = temp3;
    temp4=MODEL.A;
    MODEL.A = zeros(6*numel(MODEL.nodes(:,1)),size(temp2,2));
    MODEL.A(MODEL.freedofs,:) = temp4;
end

end

