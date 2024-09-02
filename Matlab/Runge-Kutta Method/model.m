function [dAm] = model(ts, Am, nmode, disp, dtun_ini, dtun_spd, loss_half, nlinear, Pump)
%LLE_DUALPUMPING 
%   Am(1:end/2): the envelope1 excited by pump1
%   Am(end/2+1:end): the envelope2 excited by pump2
%   ts: the slow time discrete sequence
%   disp: the dispersion coefficient of envelopei
%   pump: pump term
%   detuning_ini: initial detuning between pumpi and corresponding
%       resonance frequency
%   detuning_spd: frequency detuning speed of pumpi
%   losst: the total loss of envelopei
%   nlinear: nonlinear coefficient

%----Separate Am1, Am2 and T
    Loss_Dtun_Disper = -( loss_half + 1i*( dtun_ini + dtun_spd*ts + disp ) ).* Am ;
    At = ifft(Am).*nmode ;
    Pt = abs(At).^2 ;
    Nonlinear = 1i*nlinear.*fft(Pt.*At)./nmode ;
    dAm = Loss_Dtun_Disper + Nonlinear + Pump ;
end
