clear
clc
close all
format long
% load('GaN_wd1700_ht700_r40um.mat','p_f1','p_V1');
% load('GaN_wd1700_ht700_r40um.mat','p_f2','p_V2');
load('LN_80_1.3_air.mat','p_f','p_V');
tic;

%================================Basic Constant=======================================
c = 299792458 ;
h = 6.62e-34 ;

%================================Simulation setup=====================================
mode_span = 2^8 ;
nmode = mode_span + 1 ;
azimodes = (-mode_span/2: 1: mode_span/2)' ;
phases = 2*pi/mode_span.*(azimodes) ;
resampling_span = 4000 ;
nresampling = resampling_span + 1 ;
rephases = 2*pi/resampling_span.*(-resampling_span/2: 1: resampling_span/2)' ;

%================================Microcavity parameters===============================
Q = 1e6 ;
lamexp = 1.55e-6 ;
freqexp = c/lamexp ;
azimodeact = round(fzero(@(num) polyval(p_f,num)-freqexp, 500)) ;
freqact = polyval(p_f, azimodeact) ;
lamact = c/freqact ;
wqct = 2*pi*freqact ;
loss_t = 2*pi*freqact/Q ;
c_pct = 0.5 ;
loss_i = loss_t*(1-c_pct) ;
loss_c = loss_t*c_pct ;
ref_index = Sellmeier_Eq(lamact, '_sio2') ;
ref_index = 2.297 ;
freqs = polyval(p_f, azimodeact + azimodes) ;
wavelengths = c./freqs ;
%%%---Dispersion
pd1f = polyder(p_f) ;
pd2f = polyder(pd1f) ;
pd3f = polyder(pd2f) ;
FSR = polyval(pd1f, azimodeact) ;
trt = 1/FSR ;
D1 = 2*pi*FSR ;
D2 = 2*pi*polyval(pd2f, azimodeact) ;
D3 = 2*pi*polyval(pd3f, azimodeact) ;
Ephoton = h*freqact ;
%================================Effect in LLEs=======================================
%%----Pump
pump = 0.15 ;
Ap = sqrt(loss_c*pump/(h*freqact) ) ;
%%----Dispersion
disp = fftshift(0.5*D2*azimodes.^2 - 0*D3*azimodes.^3) ;
%%----Frequency Modulation
photon_lifetime = Q/freqact ;    % 2*pi/loss_t ;
span = 36 ;    % 15 
detuning_perphotonllife = 0.1*loss_t ;
detuning_ini = -0.1*loss_t ;
detuning_spd = detuning_perphotonllife/photon_lifetime ;
nnum = round(span*photon_lifetime/trt) ;
ts = trt.*(0:1:nnum-1);
%%----Nonlinear Term
nlinear_index = 1.8e-19 ;
Veff = polyval(p_V, azimodeact) ;
nlinear = 2*pi*h*freqact^2*nlinear_index*c/(ref_index^2*Veff) ;
%%%----Initial Field
Am = ones(nmode, 1)*exp(1i) ;

%================================Data Storage======================================================
Nstep =10;
isave = 1;
Ts = (0:1:nnum-1)*trt ;  
Ams(:, isave) = Am ;
Ats(:, isave) = ifft(Am).*nmode ;
Pstep = 500 ;
figure(1) 
plot(abs(Ats(:, isave).^2)) ;
ylim([0 5]*1e10);
Pump = Ap.*ones(nmode, 1) ;
Pumpm = zeros(nmode, 1);
Pumpm(1) = Ap ;
%================================Numerical Calculation Setup=======================================
for n=2:1:nnum
    %%----Linear Process
    L = exp( -(1i.*(disp + detuning_ini + detuning_spd*ts(n)) + 0.5*loss_t) * trt ) ;
    Am = Am.*L ;

    %%----Nonlinear Process
    At = ifft(Am).*nmode ;
    P = abs(At).^2 ;
    N = exp( 1i*nlinear.*P*trt ) ;
    At = At.*N ;
    %%----Coupling Process
    %Pump = Ap.*exp(1i*( (detuning_ini + detuning_spd*ts(n))/(2*pi) ) * phases ) ;
    Atout = 1 ;
    At = At + Pump*trt ;
    Am = fft(At)./nmode ;
%%----Monitor
    if (rem(n,Pstep) == 0)
        subplot(2, 1, 1) ;
        At = ifft(Am).*nmode ;
        plot( abs(At).^2 ) ;
        ylim([0 5]*1e10) ;
        subplot(2, 1, 2) ;
        Amsout = fftshift(Am) ;
        Pmsout = 10*log10(abs(Amsout).^2.*Ephoton*1e3) ;
        stem(freqs, Pmsout, 'Marker', 'none', 'BaseValue', -100); 
        ylim([-100 0])
        pause(0.0001) ;
    end
    %%----Save data
    if (rem(n,Nstep) == 0)
        isave = isave + 1 ;
        Ams(:, isave) = Am ;
        Ats(:, isave) = At ;
    end
    fprintf('Current progress %05.1f %% \n', ts(n)/ts(end)*100);
end
toc;
%%----Cavity state
Amssym = fftshift(Ams, 1) ;
Pmssym = Ephoton*abs(Amssym).^2 ;
PmssymdB = 10*log10(Pmssym*1e3) ;

Pts = Ephoton*abs(Ats).^2 ;

%%----Output signal
Amsout = sqrt(loss_c).*Ams ;
Amsoutsym = fftshift(Amsout, 1) ;
Pmsoutsym = Ephoton*abs(Amsoutsym).^2 ;

%================================Result Demonstration==============================================
%%----Cavity Power
figure(2);
set(gcf,'unit','normalized','position',[0.0,0.55,0.4,0.35]);
Ps = sum(Pmssym, 1) ;
subplot(2,1,1) ;
plot(Ps) ;
xlabel('Slow Time');ylabel('Power(W)');title('Cavity Power');
%%----Temperature Difference
subplot(2,1,2) ;
plot(Ts) ;
xlabel('Slow Time');ylabel('Temperature(K)');title('Temperature Difference');

%%----Specturm Evolution
figure(3) ;
set(gcf,'unit','normalized','position',[0.0,0.05,0.4,0.35]);
pcolor(PmssymdB) ;
shading flat ;

%%----Time Signal Evolution
figure(4) ;
set(gcf,'unit','normalized','position',[0.6,0.55,0.4,0.35]);
pcolor(Pts) ;
shading flat ;









