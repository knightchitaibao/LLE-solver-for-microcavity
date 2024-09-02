%****************************************************************************************************
% Normalizated field
% RKode solver
% Two envelopes' overlapping is included
% Thermal effect is included

% Riyao Zhang  2024.09.01
%****************************************************************************************************
clear
clc
close all
format long
% load('GaN_wd1700_ht700_r40um.mat','p_f1','p_V1');
% load('GaN_wd1700_ht700_r40um.mat','p_f2','p_V2');
load('LN_80_1.3_air.mat','p_f','p_V');
tic ;
%================================Basic Constant=======================================
c = 299792458 ;
h = 6.62e-34 ;

%================================Simulation setup=====================================
nmode = 2^9 ;
azimodes = (-nmode/2: nmode/2-1).' ;
phases = 2*pi/nmode.*(azimodes) ;
nresampling = 2^10 ;
rephases = 2*pi/nresampling.*(-nresampling/2: 1: nresampling/2-1)' ;
%================================Microcavity parameters===============================
%%----Cavity
r = 150e-6 ;
%%----Envelope
Q = 1e6 ;
lamexp = 1.55e-6 ;
freqexp = c/lamexp ;
azimode = round(fzero(@(num) polyval(p_f, num)-freqexp, 500)) ;
freq = polyval(p_f, azimode) ;
lam = c/freq ;
loss_t = 2*pi*freq/Q ;
loss_half = 0.5*loss_t ;
c_pct = 0.5 ;
loss_i = loss_t*(1-c_pct) ;
loss_c = loss_t*c_pct ;
ref_index = Sellmeier_Eq(lam, '_sio2') ;
ref_index = 2.297 ;
freqs = polyval(p_f, azimode + azimodes) ;
lams = c./freqs ;
lams = fftshift(lams) ;
Ephoton = h*freq ;
%%%---Dispersion
pd1f = polyder(p_f) ;
pd2f = polyder(pd1f) ;
pd3f = polyder(pd2f) ;
FSR = polyval(pd1f, azimode) ;
trt = 1/FSR ;
D1 = 2*pi*FSR ;
D2 = 2*pi*polyval(pd2f, azimode) ;
D3 = 2*pi*polyval(pd3f, azimode) ;

%================================Effect in LLEs=======================================
%%----Initial Field
Am = ones(nmode, 1)*exp(1i) ;
%%----Pump
pump = 0.16 ;
Ap = sqrt(loss_c*pump/Ephoton ) ;
Pumpm = zeros(nmode, 1) ;
Pumpm(1) = Ap ;
%%----Dispersion

%%----Frequency Modulation
photon_lifetime = Q/freq ;    % 2*pi/loss_t ;
span = 35 ;
dtun_perphotonllife = 0.15*loss_t ;
dtun_ini = -1*loss_t ;
dtun_ini_lam = dtun_ini/1.3e20 ;

dtun_spd = dtun_perphotonllife/photon_lifetime ;
dtun_tot = dtun_spd*span*photon_lifetime;
dtun_tot_lam = dtun_tot/1.3e20 ;

%%----Slow Time
ntssamp = 1000 ;
tss = photon_lifetime.*linspace(0, span, ntssamp) ;
tsendpoint = span*photon_lifetime ;
%%----Nonlinear Term
nlinear_index = 1.8e-19 ;
tsend = tss(end) ;
%================================Parameter Scanning=======================================
fin = 0.01 ;
step = 1e-4 ;
D2_ini = 2*pi*polyval(pd2f, azimode) ;
D2_fin = D2_ini*(1+fin) ;
D2_step = D2_ini*step ;
D2 = 2*pi*polyval(pd2f, azimode);
Ve_ini = polyval(p_V, azimode)*0.9 ;
Ve_fin = Ve_ini*(1+fin) ;
Ve_step = Ve_ini*step ;
options = odeset('NormControl','off',...
  'OutputFcn', @(slow_times, Am, flag) rk_monitor( slow_times, Am, flag,...
                                            nmode,...
                                            loss_c,...
                                            lams, tsend, Ephoton)) ;
m = 1 ;
%================================Numerical Calculation Setup=======================================
%%----solver setting:simultaneous monitor

clc ;
tic ;
Ve = polyval(p_V, azimode);
nlinear = 2*pi*Ephoton*freq*nlinear_index*c/(ref_index^2*Ve) ;
disper = fftshift(0.5*D2*azimodes.^2 + 0*D3*azimodes.^3) ;
[~,Ams] = ode45(@(ts, Am) model(ts, Am, nmode, disper, dtun_ini, dtun_spd, loss_half, nlinear, Pumpm),...
                                  tss, Am, options) ;
Ams = Ams.' ;
%csv_fp = sprintf('D:/Academic_Project/Dataset/Am_v1/No%d_%.4e_%.4e.csv', m, D2, Ve) ;
%csvwrite(csv_fp, Ams) ;
m = m+1 ;
toc ;
disp(m/100) ;




