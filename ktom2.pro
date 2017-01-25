;+
; NAME:
;   KTOM2
;
; PURPOSE: 
;   Determine the mass of the secondary given the velocity
;   semi-amplitude, eccentricity, inclination, Period, and mass of
;   the primary.
; 
; DESCRIPTION:
;   Determines the mass of the secondary by analytically solving the
;   cubic equation (only for the real root):
;
;   m2^3 - (m1+m2)^2*(K*sqrt(1d0-e^2)/sin(i))^3*Period/(2d0*!dpi*G) = 0
;
;   Inputs may be vectors or scalars, but must have the same number of
;   elements.
;
;   K must not be 0. 
;
; CALLING SEQUENCE:
;   m2 = ktom2(K, e, i, P, m1)
;
; INPUTS:
;   K           - The velocity semi-amplitude (m/s)
;   e           - The eccentricity of the orbit
;   i           - inclination of orbit (radians). To calculate the
;                 minimum mass (msini, assuming M1 >> M2), use
;                 i=!dpi/2d0.
;   Period      - Period of orbit (days)
;   M1          - The mass of the Primary (m_sun)
;
; OUTPUTS:
;   M2          - The mass of the secondary (m_sun)
;
; MODIFICATION HISTORY 
;  2011/08 -- Jason Eastman (Ohio State University)
;-

function ktom2, k, e, i, P, m1

G = 1.3271244d20 ;; m^3/(M_sun*second^2) Standish 1995
cubert2 = 1.25992104989487319d0
period = P*86400d0 ;; seconds

x = period/(2d0*!dpi*G)*(K*sqrt(1d0-e^2)/sin(i))^3d0
x2 = x*x
x3 = x2*x
m12 = m1*m1
m13 = m12*m1
m14 = m12*m12
;; analytic solution to cubic for the first (real) root
y = (27d0*m12*x+sqrt(729d0*m14*x2+108d0*m13*x3)+18d*m1*x2+2d*x3)^(1d/3d)
m2 = y/(3d*cubert2) - cubert2*(-6d*m1*x-x2)/(3d*y) + x/3d

return, m2

end
