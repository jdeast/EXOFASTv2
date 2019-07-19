;+
; NAME:
;   EXOFAST_GETPHASE
;
; PURPOSE:
;   Calculates the phase (mean anomaly/(2pi)) corresponding to a
;   particular true anomaly.
;
; CALLING SEQUENCE:
;   phase = exofast_getphase(eccen,omega,trueanom [,/PRIMARY]$
;                    [,/SECONDARY][,/L4][,/L5][,/PERIASTRON])
;
; EXAMPLE:
;   ;; find the phase of a primary transit described by e and omega 
;   phase = exofast_getphase(e,omega,/primary)
;
; INPUTS: 
;   ECCEN: The eccentricity of the orbit, must be between 0 and 1
;
;   OMEGA: The argument of periastron of the star, in radians
;
;   TRUEANOM: The true anomaly (in radians) for which the mean anomaly
;   is desired.
;
; OPTIONAL KEYWORDS:
;   PERIASTRON: Find the phase of Periastron (0 by definition)
;   PRIMARY: Find the phase of the primary transit
;   SECONDARY: Find the phase of the secondary eclipse
;   L4: Find the phase of the L4 Langrangian point (leading)
;   L5: Find the phase of the L5 Langrangian point (trailing)
;   ASCENDINGNODE: Find the phase of Ascending node (maximum RV signal)
;   DESCENDINGNODE: Find the phase of the Descending node (minimum RV signal)
;
; DEPENDENCIES:
;   EXOFAST_KEPLEREQ.PRO
;
; Modification History:
;  2010/06 - Rewritten by Jason Eastman (OSU)
;-

function exofast_getphase, eccen, omega, trueanom=trueanom, $
                           primary=primary, secondary=secondary, $
                           l4=l4, l5=l5, periastron=periastron,$
                           ascendingnode=ascendingnode,$
                           descendingnode=descendingnode
                           
; common phases one might be interested in
if keyword_set(periastron)     then trueanom =           0d0        ; Periastron 
if keyword_set(l5)             then trueanom =  5d0*!dpi/6d0 - omega; L5 Point (trailing)
if keyword_set(l4)             then trueanom =      !dpi/6d0 - omega; L4 Point (leading)
if keyword_set(secondary)      then trueanom =  3d0*!dpi/2d0 - omega; Secondary eclipse
if keyword_set(primary)        then trueanom =      !dpi/2d0 - omega; Primary Transit
if keyword_set(ascendingnode)  then trueanom =               - omega; Ascending Node of primary (max RV)
if keyword_set(descendingnode) then trueanom =          !dpi - omega; Descending Node of primary (min RV)

eccenanom = 2.d0*atan(sqrt((1d0-eccen)/(1d0+eccen))*tan((trueanom)/2d0))
M = eccenanom - eccen*sin(eccenanom)

phase = M/(2d0*!dpi)

neg = where(phase lt 0d0)
if neg[0] ne -1 then phase[neg] += 1d0

return, phase

end
 
