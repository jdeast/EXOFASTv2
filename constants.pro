function constants(G=G, RSun=RSun, GMSun=GMSun, day=day, GSI=GSI, c=c)

;; 2014 CODATA 
;; http://arxiv.org/abs/1507.07956
GSI = 6.67408d-11      ;; Newton's constant in SI units: m^3/(kg*s^2) 
c = 299792458d0        ;; speed of light in m/s
sigmabSI = 5.670367d-8 ;; Boltzmann's constant in J/K

;; IAU Resolution B3, 2015 
;; http://arxiv.org/abs/1510.07674
rsun = 6.957d8         ;; Radius of sun in meters
ssun = 1361d0          ;; flux at earth in W/m^2
lsun = 3.828d26        ;; luminosity of sun in W
Teffsun = 5772         ;; Effective temperature of Sun, in K
GMsun = 1.3271244d20   ;; G*Mass of sun in m^3/s^2

day = 86400d0          ;; day in seconds
rjup = 7.1492d7        ;; equatorial radius of Jupiter in m
rearth = 6.63781e6     ;; equatorial radius of Earth in m
GMearth = 3.986004d14  ;; G*Mass of Earth in m^3/s^2
GMjup = 1.2668653d17   ;; G*Mass of Jupiter in m^3/s^2

;; derived quantities in astronomically relevant units
G = GMsun/rsun^3*day^2         ;; Newton's constant in R_sun^3/(m_sun*day^2)
sigmab = sigmabSI/lsun*rsun^2  ;; Boltzmann's constant in L_sun/r_sun^2*K^4)

AU = sqrt(lsun/(4d0*!pi*ssun)) ;; AU in meters
;; Resolution B2 of the XXVIII General Assembly 
AUSI = 149597870700d0          ;; http://ssd.jpl.nasa.gov/?constants

msun = GMsun/GSI               ;; mass of sun in kg
mearth = GMearth/GSI           ;; mass of Earth in kg
mjup = GMjup/GSI               ;; mass of Earth in kg
