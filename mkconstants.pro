;+
; NAME:
;   mkconstants
; PURPOSE:
;   Returns a structure containing all the physical constants used
;   throughout EXOFASTv2, in cgs units and double precision.
;
; DESCRIPTION:
;   Unit constants (converted to cgs)
;   As defined by IAU resolutions B2, B3
;   https://www.iau.org/static/resolutions/IAU2012_English.pdf
;   https://arxiv.org/abs/1510.07674, 
;   https://arxiv.org/abs/1507.07956, Table 1
;-
function mkconstants

constants = create_struct('G', 6.67408d-8,$           ;; Gravitational Constant (cm^3/g/s^2)
                          'c', 2.99792458d10,$        ;; speed of light in vacuum (cm/s)
                          'sigmab', 5.670367d-5,$     ;; steffan boltzman constant (erg/cm^2/K^4)
                          'RSun', 6.957d10,$          ;; Solar radius (cm)
                          'LSun',3.828d33,$           ;; Solar Luminosity (erg)
                          'GMsun', 1.3271244d26,$     ;; cm^3/s^2
                          'GMearth', 3.986004d20,$    ;; cm^3/s^2
                          'GMjupiter', 1.2668653d23,$ ;; cm^3/s^2
                          'REarth',6.3781d8,$         ;; Earth Equatorial radius (cm)
                          'RJupiter',7.1492d9,$       ;; Jupiter equatorial radius (cm)
;                          'RpEarth',6.3568d8,$        ;; Earth Polar radius (cm) (not used)
;                          'RpJupiter',6.6854d9,$      ;; Jupiter Polar radius (cm) (not used)
                          'AU',1.495978707d13,$       ;; AU (cm)
                          'day',86400d0,$             ;; day (seconds)
                          'meter',100d0,$             ;; meter (cm)
                          ;; derived constants  
                          'pc',0d0,$                  ;; parsec (cm)
                          'GravitySun',0d0,$          ;; Surface gravity of sun (cm/s^2)
                          'RhoSun',0d0)               ;; Solar Density (g/cm^3)

constants.RhoSun = 3d0*constants.GMsun/constants.G/constants.RSun^3/4d0/!dpi ;; Solar Density (g/cm^3)
constants.gravitySun = constants.GMsun/constants.Rsun^2 ;; surface gravity of the sun (cm/s^2)
constants.pc = constants.au*3600d0*180d0/!dpi ;; cm in a parsec

return, constants

end
