function [Pz, normf, N_ph, Px, Py, rhoeff_normalized] = observabled2rwa(y,MM,Ns,Nb,sigmaz)
%OBSERVABLED2RWA  Norm-corrected system and bath observables for mD2.
% Extra outputs Px, Py and rhoeff_normalized are backward compatible with
% callers requesting only the historical first three outputs.

AA=reshape(y(1:MM*Ns),MM,Ns);
FF=reshape(y((MM*Ns+1):(MM*Ns+Nb*MM)),Nb,MM).';

tmpff=conj(FF)*FF.';
mff=exp(tmpff-0.5*(diag(tmpff)+diag(tmpff).'));
rhoeff = AA'*mff*AA;

normf = real(trace(rhoeff));
if ~(isfinite(normf) && normf > realmin)
    error('observabled2rwa:BadNorm', 'mD2 norm must be positive and finite.');
end
% All system observables must be divided by the variational-state norm.  The
% previous code normalized N_ph but not Pz, which could push (Pz+1)/2 outside
% [0,1] when the ODE norm drifted.
Pz = real(trace(sigmaz * rhoeff) / normf);

rhoeff_normalized = rhoeff/normf;
if Ns == 2
    sigmax = [0,1;1,0];
    sigmay = [0,-1i;1i,0];
    Px = real(trace(sigmax*rhoeff_normalized));
    Py = real(trace(sigmay*rhoeff_normalized));
else
    Px = NaN; Py = NaN;
end

% phonon in total
W = (conj(AA)*AA.').*mff;
% real part/ normf
numbath_k = real(diag(FF'*W*FF).') / normf;
N_ph = sum(numbath_k); % get sum → phonon number

end
