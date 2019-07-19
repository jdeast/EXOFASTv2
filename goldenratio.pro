;+
; NAME:
;   GOLDENRATIO
;
; PURPOSE: 
;   Uses the golden section search to find the minimum of a unimodal
;   function between the supplied bounds.  
;
;   Implementation translated from here:
;   https://en.wikipedia.org/wiki/Golden-section_search
;
; CALLING SEQUENCE:
;    min = goldenratio('func', min, max)
;
; INPUTS:
;
;    FUNC - A string specifying the name of the function to be
;           minimized
;    A    - The minimum bound on the function to be minimized
;    B    - The maximum bound on the function to be minimized
;
; OPTIONAL INPUTS:
; 
;    TOL  - The tolerance, after which the function is deemed
;           sufficiently close to the minimum
; RESULT:
;
;    The minimum value
;
; MODIFICATION HISTORY
; 
;  2018/11/13 -- Written by Jason Eastman, CfA
;
;-
function goldenratio, func, a, b, tol, use=use

if n_elements(use) eq 0 then use = lindgen(n_elements(a))
if n_elements(tol) eq 0 then tol = 1d-9

gr = (sqrt(5d0) + 1d0)/2d0

niter = 0
repeat begin

   c = b - (b - a) / gr
   d = a + (b - a) / gr

   dbig = where(call_function(func,c) lt call_function(func,d), complement=cbig)
   if dbig[0] ne -1 then b[dbig] = d[dbig]
   if cbig[0] ne -1 then a[cbig] = c[cbig]

;   print, call_function(func,c), call_function(func,d), a,b,c,d
   niter++
endrep until max(abs(c[use]-d[use])) le tol

return, (b+a)/2d0

end

