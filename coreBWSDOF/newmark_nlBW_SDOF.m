function [ MODEL ] = newmark_nlBW_SDOF( MODEL )
%Nonlinear Newmark algorithm for time integration
%
%
%Inputs:
%  MODEL : struct / The MODEL struct contains all system properties, parameters and matrices
%
%
%Returns:
%  MODEL : struct / The MODEL struct contains all system properties, parameters and matrices
%			After integrating your model forward in time the MODEL struct contains the time histories 
%			of the displacement, velocities and accelerations for every degree of freedom
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% â€œA local basis approximation approach for nonlinearparametric model order reduction,
% â€?Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

dt = MODEL.dyn.dt;
MODEL.nt=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 1/2;
beta = 1/4;

a1 = 1/(beta*dt^2);
a2 = 1/(beta*dt);
a3 = (1-2*beta)/(2*beta);
a4 = gamma/(beta*dt);
a5 = 1-gamma/beta;
a6 = (1-gamma/(2*beta))*dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = MODEL.M;
C = MODEL.dyn.a*MODEL.M + MODEL.dyn.b*MODEL.K;

a1M = a1*M; a4C = a4*C;

Fmatrix = MODEL.Rmatrix;
fnorm=norm(Fmatrix(:,MODEL.pos));
f0=Fmatrix(:,1);

v0 = MODEL.V(:,1);
if sum(MODEL.A(:,1))==0
    a0 = M\(f0-MODEL.K*MODEL.U(:,1)-C*v0);
    MODEL.A(:,1) = a0;
else
   a0= MODEL.A(:,1);
end

uk = MODEL.U(:,1);
uik = uk;

vk = v0; ak = a0; 
tol=1e-3; maxit = 20;
for i=1:(size(Fmatrix,2))
    
    fi = Fmatrix(:,i);
        
    va1 = a2*vk+a3*ak;
    va2 = a5*vk+a6*ak;
    
    R = M*(-va1) + C*va2 + MODEL.fint - fi;
    
    Rnorm = norm(R);
    
    nit=0;
    MODEL.nt = i+1;
    ui=uk;
    uik=ui-uk;
    
    while ((Rnorm>tol*fnorm)&&(nit<=maxit))
        Keff = a1M+a4C+MODEL.K;
        
        du = Keff\(-R);
        
        ui = ui+du;
        
        MODEL.u = ui;
        
        uik=ui-uk;
                
        [ MODEL ] = assemble_nlBW_SDOF( MODEL );
        [ MODEL ] = apply_bc_nl( MODEL );
        
        R = M*(a1*uik-va1) + C*(a4*uik+va2) + MODEL.fint - fi;
        
        Rnorm = norm(R);
        
        nit = nit+1;
    end
    
    if (Rnorm<tol*fnorm) 
        if (mod(i,500)==0)
            fprintf('Step %i converged after %i iterations, residual %f \n', i,nit, Rnorm/fnorm);
        end
    else
        fprintf('Step %i did not converge after %i iterations, residual %f \n', i,nit, Rnorm/fnorm);
    end
    
    MODEL.U(:,i+1)=ui;
    
    uk = ui;
    ak = a1*uik-va1;
    vk = a4*uik+va2;

    MODEL.V(:,i+1) = vk;
    MODEL.A(:,i+1) = ak;
end

end