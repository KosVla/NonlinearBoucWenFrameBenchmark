function [MODEL, Input]=LowPass(MODEL,upsamp, Input,P, N)
    if upsamp>1
        MODEL.Uups = MODEL.U;
        MODEL.U = Downsampling(MODEL.U,upsamp,P+1,N);
        MODEL.Vups = MODEL.V; MODEL.Aups = MODEL.A;
        MODEL.A = Downsampling(MODEL.A,upsamp,P+1,N);
        MODEL.V = Downsampling(MODEL.V,upsamp,P+1,N);
        f = downsample(Input.SynthesizedAccelerogram,upsamp);
        % Removal of the last simulated period to eliminate the edge effects
        %due to the low-pass filter.
        f = f(1:(P+1-1)*N,:);
        Input.SynthesizedAccelerogramUps=f;
        Input.SynthesizedAccelerogram=f;
        % P = P+1-1;
        MODEL.dyn.nt=1/f;
    end
end