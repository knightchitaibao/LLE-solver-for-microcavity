function [ri] = Sellmeier_Eq(lambda,material)
format long
lambda = lambda*1e6; % [um]
switch material
    case '_sio2'
        % McCray Handbook of Optics, Third Edition, Volume IV, Oct 2009
        ri = sqrt( 1 + 0.6961663*lambda.^2./(lambda.^2-0.0684043^2) + 0.4079426*lambda.^2./(lambda.^2-0.1162414^2) + 0.8974794*lambda.^2./(lambda.^2-9.896161^2) );%1.44624;%
    case '_sio2meas'
        % Data from CMI nk file for undensified LTO. Strangely data fits with a
        % single Sellmeier peak. No evidence for infrared absorption peak found
        % up to 1700 nm, data and fit probably not valid for wavelength larger
        % 1.8um, triple peak(fused silica) has ZDP at 1.5um but deviates from
        % data points, also infrared absorption is shifted to 223um???
        ri = sqrt( 1 + 1.0903986467543885*lambda.^2./(lambda.^2-0.091729431146042259^2) );
    case '_730'
        % LPCVD Si3N4 deposition at 730 °C @ Wu2003
        ri = sqrt( 1 + (101*0.142^2*lambda.^2)./(lambda.^2-0.142^2) );
    case '_760'
        % LPCVD Si3N4 deposition at 760 °C
        ri = sqrt( 1 + (119*0.143^2*lambda.^2)./(lambda.^2-0.143^2) );
    case '_825'
        % LPCVD Si3N4 deposition at 825°C
        ri = sqrt( 1 + (133*0.142^2*lambda.^2)./(lambda.^2-0.142^2) );
    case '_si3n4'
        % Indices from
        % http://refractiveindelambda.info/?group=CRYSTALS&material=Si3N4
        ri = sqrt( 1 + (2.8939*lambda.^2)./(lambda.^2 - 139.67e-3^2) );
    case '_si3n4geom'
        ri = 1.99;
    case '_si3n4meas'
        % Indices from spectro ellipsometric measurement fitted from 400 to
        % 1700 nm in measuredSi3N4.cfit with single peak Sellmeier Eqn.
        ri = sqrt( 1 + 2.9481912375268458*lambda.^2./(lambda.^2-0.13102873492590297^2) );
    case '_lif'
        ri = sqrt( 1 + (0.92549*lambda.^2)./(lambda.^2-0.07376^2) + (6.96747*lambda.^2)./(lambda.^2-32.79^2));
    case '_gano'
        % indicses for gallium nitride ordinary polarization from
        % http://refractiveindelambda.info/?group=CRYSTALS&material=GaN
        ri = sqrt( 3.6 + (1.75*lambda.^2)./(lambda.^2 - 0.256^2) + (4.1*lambda.^2)./(lambda.^2 - 17.86^2)  );
    case '_gane'
        % indices for gallium nitride elambdatraordinary polarization
        ri = sqrt( 5.35 + (5.08*lambda.^2)./(lambda.^2 - 17.86^2) );
    case '_hfo2'
        % from Mcray Handbook of Ooptics Vol 4 
        ri = sqrt(  1 + 1.9558*lambda.^2./(lambda.^2-0.15494^2) + 1.345*lambda.^2/(lambda.^2-0.0634^2) + 10.41*lambda.^2./(lambda.^2-27.12^2) );
    case '_zro2'
        % from Mcray Handbook of Optics Vol 4
        ri = sqrt(  1 + 1.347091*lambda.^2./(lambda.^2-0.166739^2) + 2.117788*lambda.^2./(lambda.^2-0.166739^2) + 9.452943*lambda.^2./(lambda.^2-24.320570^2) );
    case '_tio2'
        % from Mcray Handbook of Optics Vol 4
        ri = sqrt( 5.913 + 0.2441*lambda.^2./(lambda.^2-0.0803) );
    case '_si'
        % from refractiveindelambda.info
        ri = sqrt(1 + (10.6684293*lambda.^2)./(lambda.^2-0.301516485^2) + (0.003043475*lambda.^2)./(lambda.^2-1.13475115^2) + (1.54133408*lambda.^2)./(lambda.^2-1104.0^2));
    case '_mgo'
        % refractiveindelambda.info
        ri = sqrt(1+ (1.111033*lambda.^2)./(lambda.^2-0.0712465^2) + (0.8460085*lambda.^2)./(lambda.^2-0.1375204^2) + (7.808527*lambda.^2)./(lambda.^2-26.89302^2));
    case '_ge'
        ri = sqrt( 9.28156 + (6.72880*lambda.^2)./(lambda.^2-0.44105) + (0.21307*lambda.^2)./(lambda.^2-3870.1));
    case '_diamond'
        ri = sqrt( 1+ (0.3306.*lambda.^2)./(lambda.^2-0.175^2) + (4.3356*lambda.^2)./(lambda.^2-0.106^2));
    otherwise
        ri = 1;
end