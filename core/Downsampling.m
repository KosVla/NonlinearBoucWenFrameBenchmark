function [U] = Downsampling(U,upsamp,P,N)
%% Low-pass filtering and downsampling in case upsampling was employed
%
%Parameters:
%   U : matrix / Matrix to be downsamples. Every row is a signal for a dof
%   upsamp: float / Upsampling factor
%   P, N: floats / Periods and number of points per period
%
%Please cite as:
% J. Noel and M. Schoukens, 
% Hysteretic benchmark with a dynamic nonlinearity,
% in Workshop on non-linear system identification benchmarks, 2016, pp. 7–14

drate = factor(upsamp);        % prime factor decomposition.
for r=1:size(U,1)
    y=U(r,:);
    for k=1:length(drate)
        y = decimate(y,drate(k),'fir');    
    end
    y = y(1:(P-1)*N);
    U(r,1:length(y))=y;
end
U(:,length(y)+1:end) = [];

end
