function [monitor] = rk_monitor( slow_times, Am, flag,...
                                            nmode,...
                                            loss_c,...
                                            lams, tsend, Ephoton )
%   slow_times: the slow time sequence already solved
%   Am(1:end/2): the envelope1 excited by pump1
%   Am(end/2+1:end): the envelope2 excited by pump2
%   nresampling: the phase point number resampled
%   azimodes: the longitudinal mode sequence
%   dp_resampling: the phase step of resmapling
%   tsend: the slow time endpoint

monitor = 0;
disp(flag) ;
if isempty(flag)
    At = ifft(Am(1:end, end)).*nmode/(2*pi) ;
    Pm =  loss_c.*Ephoton.*abs(Am(1:end, end)).^2  ;
    Pt = Ephoton.*abs(At).^2 ;
    PmdBm = 10*log10(Pm) ;
    %--------------------------------------Representation------------------------------------------------
    set(gcf,'unit','normalized','position',[0.1,0.2,0.8,0.6]);
    %%----Time domain representation
    subplot(2,1,1);
    plot(Pt);
    ylim([0 5]*1e-11);
    %%----Wavelength representation
    subplot(2,1,2);
    stem(lams, PmdBm, 'Marker','none', 'BaseValue', -80);
    ylim([-80 0]);
    pause(0.001);
    fprintf('Current progress %05.1f %% \n', slow_times(end)/tsend*100);
end

